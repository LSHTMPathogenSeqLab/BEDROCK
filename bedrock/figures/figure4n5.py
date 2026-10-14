##genes and methylation
import pyranges as pr
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

chrom_map = {
    "Pf3D7_01_v3": "1",
    "Pf3D7_02_v3": "2",
    "Pf3D7_03_v3": "3",
    "Pf3D7_04_v3": "4",
    "Pf3D7_05_v3": "5",
    "Pf3D7_06_v3": "6",
    "Pf3D7_07_v3": "7",
    "Pf3D7_08_v3": "8",
    "Pf3D7_09_v3": "9",
    "Pf3D7_10_v3": "10",
    "Pf3D7_11_v3": "11",
    "Pf3D7_12_v3": "12",
    "Pf3D7_13_v3": "13",
    "Pf3D7_14_v3": "14",
    "Pf3D7_API_v3": "API",
    "Pf_M76611": "MIT",
}

def assign_genes(df_bed, genes_pr):
    """Add a 'gene' column. The overlap join is done once per genome position
    (not once per sample row), then merged back onto every row."""
    key = ["Chromosome", "Start_chrom_pos", "End_chrom_pos"]
    pos = df_bed[key].drop_duplicates()

    mods = pr.PyRanges(
        chromosomes=pos["Chromosome"].astype(str),
        starts=pos["Start_chrom_pos"],
        ends=pos["End_chrom_pos"],
    )

    overlaps = mods.join(genes_pr)
    odf = overlaps.df

    if odf.empty:
        raise RuntimeError(
            "No overlaps between gene annotations and BED entries."
            "Please check that chromosome names match bed files and genes."
        )

    gene_col = None

    for gene_column_name in ["Name", "gene_id", "ID", "gene_name", "locus_tag", "Parent"]:
        if gene_column_name in odf.columns:
            gene_col = gene_column_name
            break
    if gene_col is None:
        raise RuntimeError(
            f"No gene identifier column found. Avaiable columns: {odf.columns.tolist()}"
        )

    odf = odf[["Chromosome", "Start", "End", gene_col]].copy()
    odf = odf.rename(columns={gene_col: "gene", "Start": key[1], "End": key[2]})
    odf["gene"] = odf["gene"].astype("category")      # one copy of each name, fast groupby
    odf["Chromosome"] = odf["Chromosome"].astype(str).astype(df_bed["Chromosome"].dtype)
    odf[key[1]] = odf[key[1]].astype(df_bed[key[1]].dtype)
    odf[key[2]] = odf[key[2]].astype(df_bed[key[2]].dtype)

    return df_bed.merge(odf, on=key, how="left")


def summarise_gene_methylation(
    df,
    mod_codes=None,
    stat="mean",
):
    if mod_codes is not None:
        df = df[df["mod_code"].isin(mod_codes)]

    agg = {"mean": "mean", "median": "median"}[stat]

    summary = (
        df.groupby(["gene", "sample_name"], observed=True)["percent_mod"]
        .agg(agg)
        .reset_index(name=f"{stat}_methylation")
    )
    # small table: plain text labels, so plots only show genes / samples that are present
    summary["gene"] = summary["gene"].astype(str)
    summary["sample_name"] = summary["sample_name"].astype(str)

    return summary

def to_heatmap_matrix(df, value_col):
    return (
        df.pivot(index="gene", columns="sample_name", values=value_col)
        .fillna(0)
    )

def select_top_500_genes(matrix, n=500):
    variances = matrix.var(axis=1)
    top_genes = variances.sort_values(ascending=False).head(n).index
    return matrix.loc[top_genes]

def plot_heatmap(
    matrix,
    outfile,          # one path, or a list of paths (e.g. .tiff and .pdf)
    scale=None,
    cmap="vlag",
    figsize=(12, 20),
    vmin=None,
    vmax=None,
    title=None,
):
    data = matrix.copy()

    if scale == "row":
        data = (data - data.mean(axis=1).values[:, None]) / (
            data.std(axis=1).values[:, None] + 1e-6
        )

    sns.clustermap(
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        metric="euclidean",
        method="complete",
        row_cluster=False,
        col_cluster=False,
        yticklabels=False,
        figsize=figsize,
    )
    outfiles = outfile if isinstance(outfile, (list, tuple)) else [outfile]
    for f in outfiles:
        plt.savefig(f, dpi=300)
    plt.close()

def remap_chromosomes(genes_pr, chrom_map):
    df = genes_pr.df.copy()

    df["Chromosome"] = df["Chromosome"].astype(str).replace(chrom_map)

    return pr.PyRanges(df)

def run_figure4(df_bed, genes_pr, outdir):
    df_bed = df_bed[df_bed["mod_score"] >= 50] #threshold for depth
    # main.py passes rows that already have genes assigned
    df = df_bed if "gene" in df_bed.columns else assign_genes(df_bed, genes_pr)
    df = df.dropna(subset=["gene"])

    configs = [
        ("4a", None, "mean", None),
        ("4b", None, "median", "row"),
        ("4c", ["a"], "mean", None),
        ("4d", ["a"], "median", "row"),
        ("4e", ["h"], "mean", None),
        ("4f", ["h"], "median", "row"),
        ("4g", ["m"], "mean", "row"),
        ("4h", ["m"], "median", "row"),
    ]

    for label, mods, stat, scale in configs:
        summary = summarise_gene_methylation(df, mods, stat)
        if summary.empty:
            print(f"[INFO] Figure {label} skipped: no {mods} calls in genes")
            continue
        mat = to_heatmap_matrix(summary, f"{stat}_methylation")

        top_500 = select_top_500_genes(mat, n=500)

        plot_heatmap(
            top_500,
            [outdir / f"figure{label}.tiff", outdir / f"figure{label}.pdf"],
            scale=scale,
            cmap="vlag",
        )

#figure 5
##average methylation per site within gene
def run_figure5(df_bed, genes_pr, outdir):
    df_bed = df_bed[df_bed["mod_score"] >= 50]
    # main.py passes rows that already have genes assigned
    df = df_bed if "gene" in df_bed.columns else assign_genes(df_bed, genes_pr)
    df = df.dropna(subset=["gene"])

    configs = [
        ("5a", None),
        ("5b", ["a"]),
        ("5c", ["m"]),
        ("5d", ["h"]),
    ]

    for label, mods in configs:
        summary = summarise_gene_methylation(df, mods, stat="mean")
        if summary.empty:
            print(f"[INFO] Figure {label} skipped: no {mods} calls in genes")
            continue
        mat = to_heatmap_matrix(summary, "mean_methylation")

        plot_heatmap(
            mat,
            outdir / f"figure{label}.pdf",
            cmap="vlag",
            vmin=0,
            vmax=100,   # percent_mod is 0-100
            figsize=(15, 25),
        )