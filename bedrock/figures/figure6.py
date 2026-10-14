##differential methylation
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

def compute_gene_methylation(mod_summary, threshold):
    df = mod_summary.copy()

    df["methylated"] = df["mod_percentage"] >= threshold

    summary = (
        df.groupby(["gene", "sample_name"], observed=True)
        .agg(
            total_sites=("mod_percentage", "size"),
            methylated_sites=("methylated", "sum"),
        )
        .reset_index()
    )

    summary["proportion_methylated"] = (
        summary["methylated_sites"] / summary["total_sites"]
    )
    # small table: plain text labels, so plots only show genes / samples that are present
    summary["gene"] = summary["gene"].astype(str)
    summary["sample_name"] = summary["sample_name"].astype(str)

    return summary

def attach_gene_positions(methylation_df, gene_positions):
    gene_positions_unique = gene_positions.drop_duplicates("gene_id")

    out = methylation_df.merge(
        gene_positions_unique,
        left_on="gene",
        right_on="gene_id",
        how="left",
    )

    out["chr"] = out["seqnames"].astype(str)
    out["pos"] = out["midpoint"]

    return out

def plot_manhattan_by_position(df, outpath, title):
    g = sns.FacetGrid(
        df,
        row="sample_name",
        col="chr",
        sharex=False,
        sharey=True,
        height=2.5,
        aspect=1.2,
    )

    g.map_dataframe(
        sns.scatterplot,
        x="pos",
        y="proportion_methylated",
        alpha=0.7,
        s=15,
    )

    g.set_titles(row_template="{row_name}", col_template="Chr {col_name}")
    g.set_axis_labels("Genomic position", "Proportion methylated")

    for ax in g.axes.flat:
        ax.set_xticks([])
        ax.grid(False, axis="x")

    g.fig.suptitle(title, y=1.02)
    g.savefig(outpath, dpi=300)
    plt.close()

def plot_gene_manhattan(df, outpath, title):
    plt.figure(figsize=(30, 10))
    sns.scatterplot(
        data=df,
        x="gene",
        y="proportion_methylated",
        hue="sample_name",
        alpha=0.7,
        s=20,
    )

    plt.xticks(rotation=90, fontsize=5)
    plt.xlabel("Gene")
    plt.ylabel("Proportion methylated")
    plt.title(title)
    plt.grid(False, axis="x")
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()

def top_genes_overall(df, prop=0.05):
    n = int(len(df["gene"].unique()) * prop)

    top = (
        df.groupby("gene")["proportion_methylated"]
        .max()
        .nlargest(n)
        .index
    )

    return df[df["gene"].isin(top)]

def top_genes_per_sample(df, prop=0.05):
    n = int(np.ceil(len(df["gene"].unique()) * prop))
    return (
        df.sort_values("proportion_methylated", ascending=False)
        .groupby("sample_name", observed=True)
        .head(n)
    )


def plot_top_genes(df, outpath, title):
    plt.figure(figsize=(15, 10))
    sns.scatterplot(
        data=df,
        x="gene",
        y="proportion_methylated",
        hue="sample_name",
        s=25,
        alpha=0.8,
    )
    plt.xticks(rotation=90, fontsize=6)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def plot_top_genes_by_sample(df, outpath, title):
    g = sns.FacetGrid(
        df,
        col="sample_name",
        col_wrap=2,
        sharex=False,
        height=4,
    )

    g.map_dataframe(
        sns.scatterplot,
        x="gene",
        y="proportion_methylated",
        alpha=0.7,
        s=20,
    )

    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=90, labelsize=6)
        ax.grid(False, axis="x")

    g.fig.suptitle(title, y=1.02)
    g.savefig(outpath)
    plt.close()

def run_figure6(mod_summary, gene_positions, outdir):
    outdir = Path(outdir)

    for threshold, label in [(10, "10_percent"), (5, "5_percent"), (1, "1_percent")]:
        methyl = compute_gene_methylation(mod_summary, threshold)
        methyl_pos = attach_gene_positions(methyl, gene_positions)

        # 6a / 6f
        plot_manhattan_by_position(
            methyl_pos,
            outdir / f"figure6a_{label}.tiff",
            f"Manhattan plot (≥{threshold}%)",
        )

        # 6b / 6g
        plot_gene_manhattan(
            methyl,
            outdir / f"figure6b_{label}.pdf",
            f"Gene-level methylation (≥{threshold}%)",
        )

        # top 5% overall
        top_overall = top_genes_overall(methyl)
        plot_top_genes(
            top_overall,
            outdir / f"figure6c_{label}.pdf",
            "Top 5% methylated genes",
        )

        # top 5% per sample
        top_sample = top_genes_per_sample(methyl)
        plot_top_genes_by_sample(
            top_sample,
            outdir / f"figure6d_{label}.pdf",
            "Top 5% genes per sample",
        )

        # export table
        top_overall.to_csv(
            outdir / f"top5_percent_genes_methylated_{label}.csv",
            index=False,
        )