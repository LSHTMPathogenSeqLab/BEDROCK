##bases composition (failed/mod/canonical/variant)
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from bedrock.utils import make_palette
from bedrock.config import BASE_CALL_PALETTE


def base_composition(df, threshold):
    df = df[df.mod_score >= threshold]

    summary = (
        df.groupby(["Chromosome", "sample_name"], observed=True)
        [["mod", "canonical_mod", "other_mod",
          "delete_mod", "fail_mod", "variant_mod", "no_call"]]
        .sum()
        .reset_index()
    )

    long = summary.melt(
        id_vars=["Chromosome", "sample_name"],
        var_name="base_call_type",
        value_name="n_sites"
    )

    long["percent"] = (
        long.groupby(["Chromosome", "sample_name"], observed=True)
            .n_sites.transform(lambda x: x / x.sum() * 100)
    )

    long["threshold"] = f">={threshold}x"
    return long

def plot_base_composition(df, outdir):
    order = [
        "canonical_mod",
        "mod",
        "variant_mod",
        "other_mod",
        "fail_mod",
        "delete_mod",
        "no_call",
    ]

    df["base_call_type"] = pd.Categorical(
        df["base_call_type"],
        categories=order,
        ordered=True
    )

    df = (df
          .groupby(["sample_name", "threshold", "base_call_type"], as_index=False)
          .agg(percent=("percent", "sum"))
         )


    samples = df["sample_name"].unique()
    thresholds = df["threshold"].unique()
    
    base_palette = make_palette(order,palette=BASE_CALL_PALETTE)

    fig, axes = plt.subplots(
        nrows=len(samples),
        ncols=1,
        figsize=(6, 4 * len(samples)),
        sharex=True
    )


    if len(samples) == 1:
        axes = [axes]

    print(df.groupby(["sample_name", "threshold"])["base_call_type"].nunique())

    g = sns.catplot(
        data=df,
        kind="bar",
        x="threshold",
        y="percent",
        hue="base_call_type",
        col="sample_name",
        dodge=False,
        order=thresholds,
        hue_order=order,
        palette=base_palette,
        height=4,
        aspect=1
    )
    
    g.set_axis_labels("Depth Threshold", "Percentage of Base Calls")
    g.set_titles("{col_name}")
    
    g.legend.set_title("Type of Base")
    g.legend.set_bbox_to_anchor((1, 1))
    g.legend.borderaxespad=0
    g.legend.set_loc("center left")

    plt.tight_layout()
    g.savefig(f"{outdir}/figure2.tiff", dpi=300)
    g.savefig(f"{outdir}/figure2.pdf", dpi=300)
    plt.close()