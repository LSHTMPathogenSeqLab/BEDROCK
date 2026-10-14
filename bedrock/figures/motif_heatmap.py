
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import seaborn as sns

from bedrock.config import SAMPLE_PALETTE, EXTENDED_PALETTE
from bedrock.methylation_summary import MOD_NAME_MAP
from bedrock.utils import make_palette

FRACTION_VALUES = ("frac_mod", "mean_site_frac")
JOIN_COLS = ["Chromosome", "Start_chrom_pos", "strand"]


def read_table(path):
    """Read a CSV or TSV (delimiter auto-detected)."""
    df = pd.read_csv(path, sep=None, engine="python")
    df.columns = df.columns.astype(str).str.strip()
    return df


def require_columns(df, cols, what):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(
            f"[ERROR] {what} is missing column(s): {', '.join(missing)}\n"
            f"        columns found: {', '.join(df.columns)}"
        )


# ---------------------------------------------------------------- motifs
def load_motifs(motif_files, mod_code="a", required=True):
    """Union of (motif, offset) across files.

    Takes `modkit motif search` TSVs and/or any CSV/TSV with motif + offset
    columns. If a file has a mod_code column, only rows for `mod_code` are kept.
    """
    mset = set()
    for path in motif_files:
        p = Path(path)
        if not p.exists() or p.stat().st_size == 0:
            raise SystemExit(f"[ERROR] motif file not found or empty: {p}")
        df = read_table(path)
        if df.empty:          # modkit found nothing for this sample
            continue
        require_columns(df, ["motif", "offset"], f"motif file {path}")
        if "mod_code" in df.columns:
            df = df[df["mod_code"].astype(str).str.strip() == str(mod_code)]
        mset.update((str(m).strip(), int(o)) for m, o in zip(df["motif"], df["offset"]))

    if not mset:
        msg = f"No motifs for mod code '{mod_code}' in: {', '.join(map(str, motif_files))}"
        if required:
            raise SystemExit(f"[ERROR] {msg}")
        print(f"[WARN] {msg}; figure 7 is skipped")
    return sorted(mset)


def motif_site_bed(ref_fasta, motif, offset, chr_map, cache_dir):
    """All occurrences of a motif in the reference (both strands), chr-mapped.

    Built once with `modkit motif bed` and cached in cache_dir.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    bed = cache_dir / f"{motif}_{offset}.bed"

    if not bed.exists() or bed.stat().st_size == 0:
        if shutil.which("modkit") is None:
            raise SystemExit(
                "[ERROR] modkit not found on PATH; it is needed to find motif sites "
                f"(missing cache: {bed})"
            )
        tmp = bed.with_suffix(".bed.tmp")
        with open(tmp, "w") as fh:
            subprocess.run(
                ["modkit", "motif", "bed", str(ref_fasta), motif, str(offset)],
                stdout=fh, check=True,
            )
        tmp.replace(bed)

    sites = pd.read_csv(
        bed, sep="\t", header=None, usecols=[0, 1, 5],
        names=JOIN_COLS, dtype={"Chromosome": str},
    )
    sites["Chromosome"] = sites["Chromosome"].str.strip().replace(chr_map)
    sites["Start_chrom_pos"] = sites["Start_chrom_pos"].astype("int64")
    return sites.drop_duplicates()


def prepare_motif_sites(motif_files, ref_fasta, chr_map, outdir, mod_code="a", required=True):
    """Load motifs and build their site BEDs. Run early so bad input fails fast.

    Returns None when no motifs were found and required is False.
    """
    motifs = load_motifs(motif_files, mod_code, required)
    if not motifs:
        return None
    print(f"[INFO] {len(motifs)} motifs for mod code '{mod_code}': "
          + ", ".join(f"{m}({o})" for m, o in motifs))
    cache = Path(outdir) / "motif_beds"
    return {(m, o): motif_site_bed(ref_fasta, m, o, chr_map, cache) for m, o in motifs}


# ---------------------------------------------------------------- summary
def _calls_for_mod(bed, mod_code, min_cov):
    # "a" or "a,GATC,1" (pileup run with --motif) -> "a"; checked per category, not per row
    col = bed["mod_code"]
    if isinstance(col.dtype, pd.CategoricalDtype):
        wanted = [c for c in col.cat.categories if str(c).split(",")[0].strip() == str(mod_code)]
        keep = col.isin(wanted)
    else:
        keep = col.astype(str).str.split(",").str[0].str.strip() == str(mod_code)
    bm = bed.loc[keep,
                 JOIN_COLS + ["mod", "canonical_mod", "other_mod"]].copy()
    for c in ("mod", "canonical_mod", "other_mod"):
        bm[c] = pd.to_numeric(bm[c], errors="coerce").fillna(0)
    bm["Nvalid"] = bm["mod"] + bm["canonical_mod"] + bm["other_mod"]
    bm = bm[bm["Nvalid"] >= max(min_cov, 1)]
    bm["Chromosome"] = bm["Chromosome"].astype(str).str.strip()
    bm["Start_chrom_pos"] = bm["Start_chrom_pos"].astype("int64")
    return bm


def motif_methylation(samples, bed_all, motif_sites, mod_code="a", min_cov=1):
    """Long table: one row per sample x motif."""
    rows = []
    for name, bed in bed_all.groupby("sample_name", observed=True, sort=False):   # one sample at a time
        s = samples[name]
        bm = _calls_for_mod(bed, mod_code, min_cov)
        if bm.empty:
            print(f"[WARN] {name}: no '{mod_code}' calls with coverage >= {min_cov}")
        for (motif, off), sites in motif_sites.items():
            m = bm.merge(sites, on=JOIN_COLS)
            nmod, nval = m["mod"].sum(), m["Nvalid"].sum()
            rows.append({
                "sample_name": name,
                "sample_type": s["type"],
                "motif": motif,
                "offset": off,
                "label": f"{motif} ({off})",
                "n_motif_sites": len(sites),
                "n_sites": len(m),
                "Nmod": nmod,
                "Nvalid": nval,
                "frac_mod": nmod / nval if nval else np.nan,
                "mean_site_frac": (m["mod"] / m["Nvalid"]).mean() if len(m) else np.nan,
            })
    long = pd.DataFrame(rows)

    zero = long.groupby("label")["n_sites"].max() == 0
    if zero.any():
        print("[WARN] These motifs matched 0 bedMethyl sites in every sample "
              "(wrong --mod-code, chr names or reference?):")
        for lab in zero[zero].index:
            print(f"  {lab}")
    return long


def sample_metadata(samples, metadata_path=None):
    """sample_type from the spreadsheet, plus any columns from an extra CSV/TSV."""
    meta = pd.DataFrame(
        {"sample_type": {n: s["type"] for n, s in samples.items()}}
    )
    meta.index.name = "sample_name"
    if metadata_path:
        extra = read_table(metadata_path)
        # key = sample_name / sample / barcode column, else the first column
        key = next((c for c in ("sample_name", "sample", "barcode") if c in extra.columns),
                   extra.columns[0])
        extra = extra.rename(columns={key: "sample_name"})
        extra["sample_name"] = extra["sample_name"].astype(str).str.strip().str.strip("'\"")
        extra = extra.drop_duplicates("sample_name").set_index("sample_name")
        extra = extra.drop(columns=[c for c in extra.columns if c in meta.columns])
        missing = sorted(set(meta.index) - set(extra.index))
        if missing:
            print(f"[WARN] No metadata row for: {', '.join(missing)} (shown as NA)")
        meta = meta.join(extra, how="left")
    return meta


# ---------------------------------------------------------------- plotting
def plot_heatmap(long, rows, cols, value, meta, out_png, *, agg="mean",
                 sort_by=None, label_by=None, annotate=None,
                 vmin=None, vmax=None, cmap="RdYlBu_r", annot=True, title=None):
    """Pivot a long table to rows x cols and draw it, with metadata strips on top.

    meta is indexed by the `cols` values; sort_by / label_by / annotate name its columns.
    Returns the matrix that was plotted.
    """
    sort_by, label_by, annotate = list(sort_by or []), list(label_by or []), list(annotate or [])
    for c in sort_by + label_by + annotate:
        if c not in meta.columns:
            raise SystemExit(
                f"[ERROR] '{c}' is not a metadata column. Options: {', '.join(meta.columns)}"
            )

    mat = long.pivot_table(index=rows, columns=cols, values=value, aggfunc=agg)
    if mat.empty:
        raise SystemExit(f"[ERROR] Nothing to plot for '{value}'")

    order = list(mat.columns)
    if sort_by:
        key = meta.reindex(order)[sort_by].astype(str)
        key["_id"] = order
        order = key.sort_values(sort_by + ["_id"])["_id"].tolist()
    mat = mat[order]

    xlabels = [
        str(c) if not label_by
        else f"{c} (" + ", ".join(str(meta.at[c, l]) if c in meta.index else "NA" for l in label_by) + ")"
        for c in mat.columns
    ]

    if value in FRACTION_VALUES:
        vmin = 0 if vmin is None else vmin
        vmax = 1 if vmax is None else vmax

    nr, nc = mat.shape
    fig = plt.figure(figsize=(max(6, 1.1 * nc + 4), 0.45 * nr + 2.5 + 0.35 * len(annotate)))
    gs = fig.add_gridspec(
        len(annotate) + 1, 2, width_ratios=[1, 0.03], wspace=0.03, hspace=0.08,
        height_ratios=[0.35] * len(annotate) + [max(1.5, 0.45 * nr)],
    )

    for i, c in enumerate(annotate):
        ax = fig.add_subplot(gs[i, 0])
        vals = meta.reindex(mat.columns)[c].fillna("NA").astype(str)
        pal = make_palette(vals, SAMPLE_PALETTE, EXTENDED_PALETTE)   # sorted categories
        cats = list(pal)
        codes = vals.map({k: j for j, k in enumerate(cats)}).to_numpy()[None, :]
        sns.heatmap(codes, cmap=ListedColormap([pal[k] for k in cats]),
                    vmin=-0.5, vmax=len(cats) - 0.5, cbar=False,
                    xticklabels=False, yticklabels=[c], linewidths=.5, ax=ax)
        ax.tick_params(axis="y", rotation=0, length=0)
        ax.legend([Patch(color=pal[k]) for k in cats], cats, title=c,
                  loc="center left", bbox_to_anchor=(1.04, 0.5), ncol=min(len(cats), 6),
                  fontsize=8, title_fontsize=8, frameon=False)

    ax = fig.add_subplot(gs[len(annotate), 0])
    cax = fig.add_subplot(gs[len(annotate), 1])
    sns.heatmap(mat, annot=annot and nr * nc <= 600, fmt=".2f", cmap=cmap,
                vmin=vmin, vmax=vmax, linewidths=.5, ax=ax, cbar_ax=cax,
                cbar_kws={"label": value}, xticklabels=xlabels)
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.suptitle(title or f"{value} per {rows} per {cols}", y=1.0)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] {nr} x {nc} heatmap -> {out_png}")
    return mat


# ---------------------------------------------------------------- entry points
def sample_type_summary(long):
    """Pool the samples that share a sample type (spreadsheet column 2) per motif.

    frac_mod = sum(Nmod) / sum(Nvalid) over every sample of that type;
    mean_site_frac and n_sites are averaged over the samples.
    """
    grouped = (
        long.groupby(["sample_type", "motif", "offset", "label"], as_index=False, sort=False)
        .agg(n_samples=("sample_name", "nunique"),
             Nmod=("Nmod", "sum"), Nvalid=("Nvalid", "sum"),
             n_sites=("n_sites", "mean"), mean_site_frac=("mean_site_frac", "mean"))
    )
    grouped["frac_mod"] = grouped["Nmod"] / grouped["Nvalid"].where(grouped["Nvalid"] > 0)
    grouped["group"] = grouped["sample_type"] + " (n=" + grouped["n_samples"].astype(str) + ")"
    return grouped


def run_motif_heatmap(samples, bed_all, motif_sites, outdir, *, mod_code="a", min_cov=1,
                      value="frac_mod", metadata=None, sort_by=None, label_by=None,
                      annotate=None):
    """BEDROCK step: summarise motif methylation and plot two heatmaps.

    motif_heatmap_<mod>_by_sample.png       one column per sample
    motif_heatmap_<mod>_by_sample_type.png  samples with the same sample type
                                            (spreadsheet column 2) pooled into one column
    plus motif_methylation_<mod>[_by_sample_type].tsv and their _matrix.tsv files.
    Columns follow the order samples / sample types first appear in the spreadsheet.
    """
    outdir = Path(outdir)
    mod_name = MOD_NAME_MAP.get(str(mod_code), str(mod_code))

    long = motif_methylation(samples, bed_all, motif_sites, mod_code, min_cov)
    long.to_csv(outdir / f"motif_methylation_{mod_name}.tsv", sep="\t", index=False)

    # spreadsheet order, as zero-padded text so the plot's sort keeps it
    type_order = list(dict.fromkeys(s["type"] for s in samples.values()))
    type_rank = {t: f"{i:04d}" for i, t in enumerate(type_order)}

    # ---- one column per sample, grouped by sample type
    meta = sample_metadata(samples, metadata)
    meta["type_rank"] = meta["sample_type"].astype(str).map(type_rank)
    mat = plot_heatmap(
        long, "label", "sample_name", value, meta,
        outdir / f"motif_heatmap_{mod_name}_by_sample.png",
        sort_by=sort_by if sort_by is not None else ["type_rank"],
        label_by=label_by if label_by is not None else ["sample_type"],
        annotate=annotate or [],
        title=f"{mod_name} {value} per motif per sample",
    )
    mat.to_csv(outdir / f"motif_methylation_{mod_name}_matrix.tsv", sep="\t")

    # ---- one column per sample type (samples pooled)
    grouped = sample_type_summary(long)
    grouped.to_csv(outdir / f"motif_methylation_{mod_name}_by_sample_type.tsv", sep="\t", index=False)
    group_meta = grouped.drop_duplicates("group").set_index("group")[["sample_type", "n_samples"]]
    group_meta["type_rank"] = group_meta["sample_type"].map(type_rank)
    mat = plot_heatmap(
        grouped, "label", "group", value, group_meta,
        outdir / f"motif_heatmap_{mod_name}_by_sample_type.png",
        sort_by=["type_rank"],
        title=f"{mod_name} {value} per motif per sample type (samples pooled)",
    )
    mat.to_csv(outdir / f"motif_methylation_{mod_name}_by_sample_type_matrix.tsv", sep="\t")
    return long


def heatmap_from_table(table_path, rows, cols, value, outdir, *, agg="mean",
                       sort_by=None, label_by=None, annotate=None, title=None):
    """Heatmap of any long-format CSV/TSV: rows x cols coloured by value.

    Other columns become metadata for sort_by / label_by / annotate.
    """
    df = read_table(table_path)
    require_columns(df, [rows, cols, value], str(table_path))
    df[rows] = df[rows].astype(str)
    df[cols] = df[cols].astype(str)
    df[value] = pd.to_numeric(df[value], errors="coerce")
    meta = df.drop(columns=[rows, value]).groupby(cols).first()

    out = Path(outdir) / f"{Path(table_path).stem}_{value}_heatmap.png"
    mat = plot_heatmap(df, rows, cols, value, meta, out, agg=agg,
                       sort_by=sort_by, label_by=label_by, annotate=annotate, title=title)
    mat.to_csv(out.with_name(out.stem + "_matrix.tsv"), sep="\t")
    return mat