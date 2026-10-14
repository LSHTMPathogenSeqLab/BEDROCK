
import argparse
import json
import sys
from pathlib import Path

from args import add_modkit_args, add_motif_args


def get_args():
    p = argparse.ArgumentParser(
        prog="BEDROCK motif heatmap",
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    inp = p.add_argument_group("BEDROCK inputs")
    inp.add_argument("--spreadsheet", help="Same sample spreadsheet as main.py")
    inp.add_argument("--chr-list", help="Same chr map as main.py")
    inp.add_argument("--ref", help="Reference FASTA used for modkit pileup")
    p.add_argument("--outdir", default="results", help="Output directory (default: results)")

    tab = p.add_argument_group("Table mode (instead of BEDROCK inputs)")
    tab.add_argument("--table", help="Any long-format CSV/TSV to plot")
    tab.add_argument("--rows", help="Column for heatmap rows")
    tab.add_argument("--cols", help="Column for heatmap columns")
    tab.add_argument("--value", help="Column to colour by")
    tab.add_argument("--agg", default="mean", help="How to combine duplicate cells (default: mean)")

    add_modkit_args(p)
    add_motif_args(p)
    return p.parse_args()


def main():
    args = get_args()

    if args.table:
        from figures.motif_heatmap import heatmap_from_table
        if not (args.rows and args.cols and args.value):
            sys.exit("[ERROR] --table needs --rows, --cols and --value")
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        heatmap_from_table(
            args.table, args.rows, args.cols, args.value, outdir,
            agg=args.agg,
            sort_by=args.motif_sort_by,
            label_by=args.motif_label_by,
            annotate=args.motif_annotate,
        )
        return

    from check_install import check_install
    from data_loading import load_chr_map, load_samples, concat_samples
    from figures.motif_heatmap import prepare_motif_sites, run_motif_heatmap
    from modkit_runner import prepare_inputs, find_motif_tables

    check_install()
    for name in ("spreadsheet", "chr_list", "ref"):
        path = getattr(args, name)
        if not path:
            sys.exit(f"[ERROR] --{name.replace('_', '-')} is required (or use --table)")
        p = Path(path)
        if not p.exists() or p.stat().st_size == 0:
            sys.exit(f"[ERROR] {name} file not found or empty: {p}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "run_args.json", "w") as fh:
        json.dump(vars(args), fh, indent=2)
    chr_map = load_chr_map(args.chr_list)

    spreadsheet = prepare_inputs(
        args.spreadsheet, args.ref, outdir,
        threads=args.threads, pileup_args=args.pileup_args, rerun=args.rerun_modkit,
    )
    motif_files = args.motifs or find_motif_tables(
        spreadsheet, args.ref, outdir, mod_code=args.mod_code, threads=args.threads,
        extra_args=args.motif_search_args, rerun=args.rerun_modkit,
    )
    if not motif_files:
        sys.exit("[ERROR] No motif tables: install modkit or give --motifs")
    motif_sites = prepare_motif_sites(motif_files, args.ref, chr_map, outdir, args.mod_code)
    samples = load_samples(spreadsheet, chr_map)
    bed_all = concat_samples(samples, "bed")

    run_motif_heatmap(
        samples,
        bed_all,
        motif_sites,
        outdir,
        mod_code=args.mod_code,
        min_cov=args.motif_min_cov,
        value=args.motif_value,
        metadata=args.sample_metadata,
        sort_by=args.motif_sort_by,
        label_by=args.motif_label_by,
        annotate=args.motif_annotate,
    )
    print("[INFO] Motif heatmap produced")


if __name__ == "__main__":
    main()