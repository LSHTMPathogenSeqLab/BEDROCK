##create figure 1 for QC = depth/coverage
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from bedrock.utils import make_palette
from bedrock.config import SAMPLE_PALETTE

def figure_1a(samples, outdir):
    df = pd.concat(
        [
            s["depth"].assign(sample=name)
            for name, s in samples.items()
        ]
    )

    sample_palette = make_palette(
        df["sample"].unique(),
        palette=SAMPLE_PALETTE
    )

    summary = (
        df.groupby(["Chromosome", "sample"], observed=True)
          .Coverage_depth.mean()
          .reset_index()
    )

    plt.figure(figsize=(8, 6))
    sns.barplot(
        data=summary,
        x="Chromosome",
        y="Coverage_depth",
        hue="sample",
        palette=sample_palette
    )
    plt.yscale("log")
    plt.ylabel("Mean Coverage Depth (log10)")
    plt.tight_layout()
    plt.savefig(f"{outdir}/figure_1a.tiff", dpi=300)
    plt.savefig(f"{outdir}/figure_1a.pdf", dpi=300)
    plt.close()


def figure_1b(df_depth_all, fai, thresholds, outdir):
    out = []

    for t in thresholds:
        tmp = (
            df_depth_all[df_depth_all.Coverage_depth >= t]
            .groupby(["Chromosome", "sample_name"], observed=True)
            .size()
            .reset_index(name="covered")
            .merge(fai, on="Chromosome")
        )
        tmp["percent"] = tmp.covered / tmp.length * 100
        tmp["threshold"] = f">={t}x"
        out.append(tmp)

    df = pd.concat(out)

    sample_palette = make_palette(
        df["sample_name"].unique(),
        palette=SAMPLE_PALETTE
    )


    g = sns.catplot(
        data=df,
        x="Chromosome",
        y="percent",
        hue="sample_name",
        col="threshold",
        kind="bar",
        col_wrap=3,
        height=4,
        palette=sample_palette
    )
    g.set_axis_labels("Chromosome", "Percent Covered")
    g.savefig(f"{outdir}/figure_1b.tiff", dpi=300)
    g.savefig(f"{outdir}/figure_1b.pdf", dpi=300)
