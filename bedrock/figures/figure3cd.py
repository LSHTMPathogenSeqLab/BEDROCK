#modified bases by AT or CG in 1kb
import pandas as pd
import numpy as np
import pyranges as pr
import seaborn as sns
import matplotlib.pyplot as plt

from bedrock.utils import make_palette

WINDOW_SIZE = 1000

def get_chromosome_sizes(bed_df):

    chrom_sizes = (
        bed_df
        .groupby("Chromosome", observed=True)["End_chrom_pos"]
        .max()
        .reset_index()
        .rename(columns={"End_chrom_pos": "chrom_length"})
    )
    
    return chrom_sizes

def make_1kb_windows(chrom_sizes, window_size = 1000):
    windows = []
    for _, row in chrom_sizes.iterrows():
        chrom = row["Chromosome"]
        length = row["chrom_length"]

        starts = np.arange(0, length, window_size)
        ends = starts + window_size
        ends[ends > length] = length

        df = pd.DataFrame({
            "Chromosome": chrom,
            "Start": starts,
            "End": ends
        })

        windows.append(df)

    windows_df = pd.concat(windows, ignore_index=True)

    return pr.PyRanges(windows_df)

def bed_to_ranges(df_bed):
    df = df_bed.rename(columns={
        "Start_chrom_pos": "Start",
        "End_chrom_pos": "End"
    })
    
    return pr.PyRanges(df)


def summarize_modifications(bed_df, window_pr, mod_codes):
    # filter first, then copy/convert only the rows that are used
    keep = (bed_df["mod_code"].isin(mod_codes)) & (bed_df["mod"] == 1)
    filtered = bed_df.loc[keep, ["Chromosome", "Start_chrom_pos", "End_chrom_pos",
                                 "strand", "sample_name"]].rename(columns={"sample_name": "sample"})
    if filtered.empty:
        print(f"[INFO] No {mod_codes} modifications found; skipping that panel")
        return pd.DataFrame()

    filtered["sample"] = filtered["sample"].astype(str).str.strip()
    filtered["strand"] = filtered["strand"].astype(str).str.strip()
    filtered["Chromosome"] = filtered["Chromosome"].astype(str).str.strip()
    
    bed = bed_to_ranges(filtered)
    overlaps = bed.join(window_pr)
    df = overlaps.df

    df["sample"] = df["sample"].astype(str).str.strip()
    df["strand"] = df["strand"].astype(str).str.strip()
    df["Chromosome"] = df["Chromosome"].astype(str).str.strip()

    summary = (
        df.groupby(
            ["Chromosome", "Start_b", "End_b", "strand", "sample"], as_index=False
        )
        .size()
        .rename(columns={
            "size": "mod_count",
            "Start_b": "window_start",
            "End_b": "window_end"
        })
    )

    strand_map = {"+": 1, "-": -1}
    summary["mod_count"] = (summary["mod_count"] * summary["strand"].map(strand_map))

    all_windows = (window_pr.df.rename(columns={"Start": "window_start", "End": "window_end"}).copy())
    all_windows["Chromosome"] = all_windows["Chromosome"].astype(str)

    samples = filtered["sample"].unique()
    strands = ["+", "-"]

    template = all_windows.merge(
            pd.MultiIndex.from_product(
                [samples, strands],
                names=["sample", "strand"]
            ).to_frame(index=False),
            how="cross"
    )

    template["sample"] = template["sample"].astype(str)
    template["strand"] = template["strand"].astype(str)
    template["Chromosome"] = template["Chromosome"].astype(str)

    summary["sample"] = summary["sample"].astype(str)
    summary["strand"] = summary["strand"].astype(str)
    summary["Chromosome"] = summary["Chromosome"].astype(str)

    merged = template.merge(
        summary,
        on=["Chromosome", "window_start", "window_end", "strand", "sample"],
        how="left"
    )
    merged["mod_count"] = merged["mod_count"].fillna(0)
    
    return merged
    
    
def plot_mod_windows(df, outpath, ylab):

    if df.empty:
        print(f"Skipping plot {outpath}:no data. Might be an error or your sample does not contain the specific modification.")
        return

    print(df.columns)
    print(df["mod_count"].abs().max())

    df = df.copy()

    chrom_order = list(dict.fromkeys(df["Chromosome"]))
    df["Chromosome"] = pd.Categorical(
        df["Chromosome"],
        categories=chrom_order,
        ordered=True
    )

    present_samples = sorted(df["sample"].dropna().unique())
    df["sample"] = pd.Categorical(
        df["sample"],
        categories=present_samples,
        ordered=True
    )

    g = sns.FacetGrid(
        df,
        row="sample",
        col="Chromosome",
        sharex=False,
        sharey=True,
        height=3
    )

    def draw(data, **kwargs):
        ax = plt.gca()

        for strand in ["+", "-"]:
            sub = data[data["strand"] == strand].sort_values("window_start")
            color = "red" if strand == "+" else "blue"
            starts = sub["window_start"].to_numpy()
            ends = sub["window_end"].to_numpy()

            if len(sub) and (starts[1:] == ends[:-1]).all():
                # contiguous windows: one filled step outline instead of one bar per window
                ax.stairs(sub["mod_count"].to_numpy(), np.append(starts, ends[-1]),
                          baseline=0, fill=True, color=color)
            else:
                ax.bar(
                    sub["window_start"],
                    sub["mod_count"],
                    width=WINDOW_SIZE,
                    align="edge",
                    color=color
                )

        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels((ax.get_xticks() / 1000).astype(int))

    g.map_dataframe(draw)
    g.set_axis_labels("Genomic Position (kb)", ylab)
    g.set_titles(col_template="{col_name}", row_template="")

    num_cols = len(df["Chromosome"].cat.categories)
    max_name_len = max(len(str(s)) for s in g.row_names)
    right_margin = max(0.85, 1 - (max_name_len * 0.008))
    
    for i, sample_name in enumerate(g.row_names):
        ax = g.axes[i, num_cols - 1]
        bbox = ax.get_position()
        g.figure.text(
            right_margin + 0.01,
            bbox.y0 + bbox.height / 2,
            str(sample_name).strip("'\""),
            va="center",
            ha="left",
            fontsize=10,
            transform=g.figure.transFigure
        )

    plt.subplots_adjust(hspace=0.4, right=right_margin)
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()