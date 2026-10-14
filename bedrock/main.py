from pathlib import Path
import pandas as pd
import sys
import json
import re
import resource
import time

from args import get_args
from data_loading import load_chr_map, load_fai, load_samples, load_genes_as_pyranges, concat_samples
from check_install import check_install
from methylation_summary import methylation_summary
from reference_counts import count_reference_bases
from figures.figure1 import figure_1a, figure_1b
from figures.figure2 import base_composition, plot_base_composition
from figures.merge_figure3 import run_figure3
from figures.figure4n5 import run_figure4, assign_genes, run_figure5, remap_chromosomes
from figures.figure6 import run_figure6
from figures.motif_heatmap import prepare_motif_sites, run_motif_heatmap
from modkit_runner import prepare_inputs, find_motif_tables

_T0 = [time.perf_counter()]


def tick(step):
    """Print how long the step since the last tick took, and peak memory so far."""
    now = time.perf_counter()
    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2   # KB on Linux
    print(f"[TIME] {step}: {now - _T0[0]:.1f}s (peak memory {peak_gb:.1f} GB)", flush=True)
    _T0[0] = now


def main():
    print(f"BEDROCK - Bedmethyl Visualisation Tool; V0.1; used Modkit output only")

    # ---- check packages ----
    check_install()
    
    args = get_args()

    # ---- validate argument inputs ----

    required_files = {
        "spreadsheet": args.spreadsheet,
        "chr_list": args.chr_list,
        "ref": args.ref,
    }

    for name, path in required_files.items():
        p = Path(path)
        if not p.exists():
            sys.exit(f"[ERROR] {name} file not found: {p}")
        if p.stat().st_size == 0:
            sys.exit(f"[ERROR] {name} file is empty: {p}")
    

    # ---- create output directory ----
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(outdir / "run_args.json", "w") as fh:
        json.dump(vars(args), fh, indent=2)


    # ---- load data ----
    chr_map = load_chr_map(args.chr_list)
    fai = load_fai(f"{args.ref}.fai", chr_map)
    thresholds = args.coverage_thresholds

    # ---- modkit: BAM rows -> bedMethyl (pileup) + depth ----
    spreadsheet = prepare_inputs(
        args.spreadsheet, args.ref, outdir,
        threads=args.threads, pileup_args=args.pileup_args, rerun=args.rerun_modkit,
    )

    # ---- motifs for figure 7 (before the long figures, so problems show early) ----
    motif_sites = None
    if not args.skip_motifs:
        motif_files = args.motifs or find_motif_tables(
            spreadsheet, args.ref, outdir, mod_code=args.mod_code, threads=args.threads,
            extra_args=args.motif_search_args, rerun=args.rerun_modkit,
        )
        if motif_files:
            motif_sites = prepare_motif_sites(
                motif_files, args.ref, chr_map, outdir, args.mod_code,
                required=bool(args.motifs),
            )

    tick("modkit / motif setup")
    samples = load_samples(spreadsheet, chr_map)
    tick("load bedMethyl + depth")

    # one table per data type for all samples; per-sample bedMethyl tables are
    # released after this, so the (large) methylation data is only held once
    df_depth_all = concat_samples(samples, "depth")
    bed_all = concat_samples(samples, "bed")
    tick("combine samples")

    # ---- FIGURE 1 ----
    figure_1a(
        samples,
        outdir
    )

    print(f"[INFO] Figure 1a produced")

    tick("Figure 1a")

    figure_1b(
        df_depth_all,
        fai,
        thresholds,
        outdir,
    )

    print(f"[INFO] Figure 1b produced")

    tick("Figure 1b")

    # ---- FIGURE 2 ----

    df_base = pd.concat(
        base_composition(bed_all, t) for t in thresholds
    )

    plot_base_composition(
        df_base,
        outdir
    )

    print(f"[INFO] Figure 2 produced")

    tick("Figure 2")

    # ---- Methylation Summary ----
    methylation_summary(
        bed_all=bed_all,
        outdir=outdir,
    )
    print(f"[INFO] Methylation summary produced")
    tick("Methylation summary")

    # ---- FIGURE 3A/B ----
    run_figure3(
        ref_fasta=args.ref,
        bed_all=bed_all,
        chr_map=chr_map,
        outdir=outdir,
    )

    print(f"[INFO] Figure 3 produced")

    tick("Figure 3")

    # ---- FIGURE 7: motif methylation heatmap (needs --motifs) ----
    if motif_sites:
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
        print(f"[INFO] Figure 7 produced")
        tick("Figure 7")
    else:
        print(f"[INFO] Figure 7 skipped (no motifs)")

    if args.annotation is None:
        print("[INFO] No GFF provided: skipping figures 4, 5 and 6. Finished.")
        return

    # ---- FIGURE 4&5 ----
    genes_pr = load_genes_as_pyranges(args.annotation)
    genes_pr = remap_chromosomes(genes_pr, chr_map)

    # assign genes once (used by figures 4, 5 and 6), keeping only the columns they use
    df_gene = assign_genes(
        bed_all[["Chromosome", "Start_chrom_pos", "End_chrom_pos", "mod_code",
                 "mod_score", "percent_mod", "sample_name"]],
        genes_pr,
    )
    df_gene = df_gene.dropna(subset=["gene"])
    del bed_all          # figures 4-6 only need the gene-assigned rows
    tick("assign genes")

    run_figure4(
        df_gene, 
        genes_pr, 
        outdir,
    )
    
    print(f"[INFO] Figure 4 produced")
    
    tick("Figure 4")

    run_figure5(
        df_gene, 
        genes_pr, 
        outdir,
    )
    
    print(f"[INFO] Figure 5 produced")
    
    tick("Figure 5")

    # ---- FIGURE 6 ----

    # clean gene names once per gene (category), not once per row
    def _clean(name):
        return re.sub(r"^gene:", "", str(name).strip().strip("'\""))

    genes_col = df_gene["gene"].astype("category")
    cleaned = [_clean(c) for c in genes_col.cat.categories]
    if len(set(cleaned)) == len(cleaned):
        df_gene["gene"] = genes_col.cat.rename_categories(cleaned)
    else:
        df_gene["gene"] = genes_col.astype(str).map(_clean).astype("category")

    mod_summary = (
        df_gene[["gene", "sample_name", "percent_mod"]]
        .rename(columns={"percent_mod": "mod_percentage"})
    )

    gene_positions = (
        genes_pr.df
        .assign(midpoint=lambda x: (x.Start + x.End) // 2,
        gene_id=lambda x: x["gene_id"].astype(str).str.replace("^gene:", "", regex=True))
        [["gene_id", "Chromosome", "midpoint"]]
        .rename(columns={"Chromosome": "seqnames"})
    )

    run_figure6(
        mod_summary,
        gene_positions,
        outdir,
    )

    print(f"[INFO] Figure 6 produced")

    tick("Figure 6")

if __name__ == "__main__":
    main()