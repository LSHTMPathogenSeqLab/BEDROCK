##AT or CG pairs per 1kb over forward and reverse strands
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO

from bedrock.utils import make_palette


WINDOW_SIZE = 1000

##count bases in FR strands - fixed window size
def window_base_counts(fasta_path, chrom_map, bases):

    records = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom = chrom_map.get(record.id, record.id)
        seq = str(record.seq).upper()
        rev = str(record.seq.reverse_complement()).upper()

        for start in range(0, len(seq), WINDOW_SIZE):
            end = min(start + WINDOW_SIZE, len(seq))

            for base, strand_seq, strand in [
                (bases[0], seq, "+"),
                (bases[1], rev, "-"),
            ]:
                count = strand_seq[start:end].count(base)
                records.append(
                    dict(
                        Chromosome=chrom,
                        start=start,
                        end=end,
                        strand=strand,
                        Count=count if strand == "+" else -count,
                    )
                )

    return pd.DataFrame(records)


def plot_window_counts(df, outpath, ylab):

    df["strand"] = pd.Categorical(
        df["strand"],
        categories=["+", "-"],
        ordered=True
    )

    chrom_order = list(dict.fromkeys(df["Chromosome"]))

    df["Chromosome"] = pd.Categorical(
        df["Chromosome"],
        categories=chrom_order,
        ordered=True
    )

    palette = make_palette(df["strand"].unique(), palette="Set1")

    g = sns.FacetGrid(
        df,
        col="Chromosome",
        col_wrap=4,
        sharex=False,
        sharey=True,
        height=3,
    )

    g.map_dataframe(
        sns.barplot,
        x="start",
        y="Count",
        hue="strand",
        palette=palette,
        dodge=False,
    )

    g.add_legend()
    g.set_axis_labels("Genomic position (bp)", ylab)

    for ax in g.axes.flatten():
        ax.tick_params(axis="x", labelbottom=False)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    
    plt.close()
