# run modkit (pileup, motif search) and samtools depth so BEDROCK can start from modBAMs
#
# Spreadsheet column 3 may be a modBAM/CRAM instead of a bedMethyl:
#   BAM  -> modkit pileup         -> <outdir>/bedmethyls/<sample>.bedmethyl
# Spreadsheet column 4 may be "-" (or empty) when column 3 is a BAM:
#   BAM  -> samtools depth -a     -> <outdir>/depth/<sample>.depth
# Figure 7 motifs, unless --motifs is given:
#   bedMethyl -> modkit motif search -> <outdir>/motifs/<sample>.motifs.tsv
#
# Every output is reused on the next run if it is newer than its input
# (--rerun-modkit forces a rerun).
import re
import shlex
import shutil
import subprocess
import sys
from functools import lru_cache
from pathlib import Path

import pandas as pd

SHEET_COLS = ["sample_name", "sample_type", "bedmethyl_path", "depth_path"]
ALN_EXT = (".bam", ".cram")
NO_FILE = {"", "-", "na", "nan", "none"}
MIN_MODKIT = (0, 4)   # `modkit motif search` / `modkit motif bed` layout


# ---------------------------------------------------------------- helpers
def require_tool(name, why):
    if shutil.which(name) is None:
        sys.exit(f"[ERROR] '{name}' not found on PATH; it is needed to {why}")


@lru_cache(maxsize=None)
def check_modkit():
    """Exit unless modkit >= 0.4 is on PATH. Returns the version string."""
    require_tool("modkit", "run modkit pileup / motif search")
    out = subprocess.run(["modkit", "--version"], capture_output=True, text=True).stdout
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", out)
    if not m:
        sys.exit(f"[ERROR] could not read the modkit version from: {out.strip()!r}")
    if (int(m.group(1)), int(m.group(2))) < MIN_MODKIT:
        sys.exit(f"[ERROR] BEDROCK needs modkit >= {'.'.join(map(str, MIN_MODKIT))}; found {m.group(0)}")
    print(f"[INFO] Using modkit {m.group(0)}")
    return m.group(0)


def _run(cmd, stdout=None):
    cmd = [str(c) for c in cmd]
    print("[CMD] " + " ".join(shlex.quote(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True, stdout=stdout)


def _up_to_date(out, src, rerun):
    out = Path(out)
    return (not rerun and out.exists() and out.stat().st_size > 0
            and out.stat().st_mtime >= Path(src).stat().st_mtime)


def _tmp(path):
    tmp = Path(path).with_name(Path(path).name + ".tmp")
    tmp.unlink(missing_ok=True)
    return tmp


def _safe(name):
    return re.sub(r"[^\w.-]+", "_", name)


def _is_alignment(path):
    return Path(path).suffix.lower() in ALN_EXT


# ---------------------------------------------------------------- tools
def ensure_index(aln):
    aln = Path(aln)
    if any(Path(str(aln) + e).exists() for e in (".bai", ".csi", ".crai")) or aln.with_suffix(".bai").exists():
        return
    require_tool("samtools", "index alignment files")
    _run(["samtools", "index", aln])


def run_pileup(aln, ref, out_bed, threads=4, extra_args="", rerun=False):
    """modBAM -> bedMethyl (all tab-separated)."""
    out_bed = Path(out_bed)
    if _up_to_date(out_bed, aln, rerun):
        print(f"[INFO] Reusing {out_bed}")
        return out_bed
    check_modkit()
    ensure_index(aln)
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp(out_bed)
    _run(["modkit", "pileup", aln, tmp, "--ref", ref, "--threads", threads,
          "--log-filepath", out_bed.with_suffix(".log")] + shlex.split(extra_args))
    tmp.replace(out_bed)
    return out_bed


def run_depth(aln, out_depth, ref=None, rerun=False):
    """samtools depth -a (every position, 3 columns)."""
    out_depth = Path(out_depth)
    if _up_to_date(out_depth, aln, rerun):
        print(f"[INFO] Reusing {out_depth}")
        return out_depth
    require_tool("samtools", "calculate depth")
    out_depth.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["samtools", "depth", "-a", aln]
    if Path(aln).suffix.lower() == ".cram" and ref:
        cmd[2:2] = ["--reference", ref]
    tmp = _tmp(out_depth)
    with open(tmp, "w") as fh:
        _run(cmd, stdout=fh)
    tmp.replace(out_depth)
    return out_depth


def run_motif_search(bedmethyl, ref, out_tsv, mod_code="a", threads=4, extra_args="", rerun=False):
    """bedMethyl -> table of modification-enriched motifs. Returns None if the search fails."""
    out_tsv = Path(out_tsv)
    if _up_to_date(out_tsv, bedmethyl, rerun):
        print(f"[INFO] Reusing {out_tsv}")
        return out_tsv
    check_modkit()
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    tmp = _tmp(out_tsv)
    try:
        _run(["modkit", "motif", "search", "-i", bedmethyl, "-r", ref, "-o", tmp,
              "--mod-code", mod_code, "--threads", threads, "--suppress-progress",
              "--log-filepath", out_tsv.with_suffix(".log")] + shlex.split(extra_args))
    except subprocess.CalledProcessError as e:
        print(f"[WARN] modkit motif search failed for {bedmethyl} (exit {e.returncode}); "
              f"see {out_tsv.with_suffix('.log')}")
        return None
    if not tmp.exists():
        print(f"[WARN] modkit motif search wrote no table for {bedmethyl}")
        return None
    tmp.replace(out_tsv)
    return out_tsv


# ---------------------------------------------------------------- pipeline steps
def read_sheet(spreadsheet):
    sheet = pd.read_csv(spreadsheet, sep="\t", header=None, names=SHEET_COLS,
                        dtype=str, keep_default_na=False, comment=None)
    sheet = sheet.fillna("")
    for c in SHEET_COLS:
        sheet[c] = sheet[c].astype(str).str.strip().str.strip("'\"")
    return sheet[sheet["sample_name"] != ""]


def prepare_inputs(spreadsheet, ref, outdir, threads=4, pileup_args="", rerun=False):
    """Turn BAM rows into bedMethyl + depth files.

    Returns the spreadsheet to load: the original one if every row already
    points at a bedMethyl and depth file, else <outdir>/samples_resolved.tsv.
    """
    outdir = Path(outdir)
    sheet = read_sheet(spreadsheet)
    changed = False
    rows = []

    for r in sheet.itertuples(index=False):
        src = Path(r.bedmethyl_path)
        if not src.exists():
            sys.exit(f"[ERROR] {r.sample_name}: file in column 3 not found: {src}")

        if _is_alignment(src):
            bed = run_pileup(src, ref, outdir / "bedmethyls" / f"{_safe(r.sample_name)}.bedmethyl",
                             threads, pileup_args, rerun)
            changed = True
        else:
            bed = src

        if r.depth_path.lower() in NO_FILE:
            if not _is_alignment(src):
                sys.exit(f"[ERROR] {r.sample_name}: no depth file in column 4. Give one, or put the "
                         "BAM in column 3 so BEDROCK can run samtools depth")
            depth = run_depth(src, outdir / "depth" / f"{_safe(r.sample_name)}.depth", ref, rerun)
            changed = True
        else:
            depth = Path(r.depth_path)
            if not depth.exists():
                sys.exit(f"[ERROR] {r.sample_name}: depth file not found: {depth}")

        rows.append([r.sample_name, r.sample_type, str(Path(bed).resolve()), str(Path(depth).resolve())])

    if not changed:
        return spreadsheet
    resolved = outdir / "samples_resolved.tsv"
    pd.DataFrame(rows).to_csv(resolved, sep="\t", header=False, index=False)
    print(f"[INFO] bedMethyl/depth paths used for this run: {resolved}")
    return resolved


def find_motif_tables(spreadsheet, ref, outdir, mod_code="a", threads=4, extra_args="", rerun=False):
    """Run modkit motif search on every sample's bedMethyl. Returns the tables that were made."""
    if shutil.which("modkit") is None:
        print("[WARN] modkit not on PATH: cannot search for motifs, so figure 7 is skipped "
              "(install modkit or give --motifs)")
        return []
    tables = []
    for r in read_sheet(spreadsheet).itertuples(index=False):
        t = run_motif_search(r.bedmethyl_path, ref,
                             Path(outdir) / "motifs" / f"{_safe(r.sample_name)}.motifs.tsv",
                             mod_code, threads, extra_args, rerun)
        if t is not None:
            tables.append(t)
    return tables
