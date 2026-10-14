import argparse


def add_modkit_args(parser):
    """Options for running modkit / samtools on BAM rows of the spreadsheet."""
    g = parser.add_argument_group(
        "modkit (used when spreadsheet column 3 is a modBAM, or to find motifs)"
    )
    g.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for modkit pileup / motif search (default: 4)",
    )
    g.add_argument(
        "--pileup-args",
        default="",
        help='Extra options for modkit pileup, in quotes, e.g. "--filter-threshold 0.7"',
    )
    g.add_argument(
        "--motif-search-args",
        default="",
        help='Extra options for modkit motif search, e.g. "--min-sites 100 --search-timeout 30m"',
    )
    g.add_argument(
        "--rerun-modkit",
        action="store_true",
        help="Rerun pileup / depth / motif search even if their outputs are already in --outdir",
    )
    return parser


def add_motif_args(parser):
    """Options for the motif heatmap (Figure 7). Shared by main.py and run_motif_heatmap.py."""
    g = parser.add_argument_group("Motif heatmap (Figure 7)")
    g.add_argument(
        "--motifs",
        nargs="+",
        metavar="FILE",
        help="Use these motif tables (modkit `motif search` TSVs or any CSV/TSV with motif + "
             "offset columns) instead of running modkit motif search on every sample",
    )
    g.add_argument(
        "--skip-motifs",
        action="store_true",
        help="Do not make figure 7 (no motif search)",
    )
    g.add_argument(
        "--mod-code",
        default="a",
        help="Modification to summarise: a=6mA, m=5mC, h=5hmC, 21839=4mC (default: a)",
    )
    g.add_argument(
        "--motif-min-cov",
        type=int,
        default=1,
        help="Ignore motif sites with fewer valid reads than this (default: 1)",
    )
    g.add_argument(
        "--motif-value",
        default="frac_mod",
        choices=["frac_mod", "mean_site_frac", "n_sites"],
        help="Value to colour the heatmap by (default: frac_mod)",
    )
    g.add_argument(
        "--sample-metadata",
        help="Optional CSV/TSV of extra sample info (lineage, country, ...). "
             "Matched on a sample_name / sample / barcode column, else the first column",
    )
    g.add_argument(
        "--motif-sort-by",
        nargs="+",
        help="Metadata column(s) to order samples by (default: sample_type)",
    )
    g.add_argument(
        "--motif-label-by",
        nargs="+",
        help="Metadata column(s) added to the sample labels",
    )
    g.add_argument(
        "--motif-annotate",
        nargs="+",
        help="Metadata column(s) drawn as coloured strips above the heatmap (default: sample_type)",
    )
    return parser


def get_args():
    parser = argparse.ArgumentParser(
        prog="BEDROCK",
        description="BEDROCK – Bedmethyl visualisation tool (Modkit output only)",
        epilog=(
            "Example usage:\n"
            "  python main.py "
            "--spreadsheet samples.tsv "
            "--chr-list chr_map.tsv "
            "--ref ref.fasta "
            "--annotation genes.gff3 "
            "--coverage-thresholds 50 100 150 200 250 300 "
            "--outdir results\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--spreadsheet",
        required=True,
        help="TSV, no header: sample_name, sample_type, bedMethyl OR modBAM path, "
             "depth path (or - to compute it from the BAM)"
    )

    parser.add_argument(
        "--chr-list",
        required=True,
        help="Chr list should have #1=ref_chr;ref_chr name #2=assigned_chr;assigned chr number/name"
    )

   

    parser.add_argument(
        "--ref",
        required=True,
        help="Reference FASTA file (.fasta)"
    )

    parser.add_argument(
        "--annotation",
        required=False,
        help="GFF3 annotation file (.gff); optional, figures 4-6 are skipped without it"
    )

    parser.add_argument(
        "--outdir",
        default="results",
        help="Output directory (default: results)"
    )

    parser.add_argument(
        "--coverage-thresholds",
        type=int,
        nargs="+",
        default=[50, 100, 150, 200, 250, 300],
        help="Coverage thresholds for QC and composition plots"
    )

    parser.add_argument(
        "--no-log-depth",
        action="store_true",
        help="Disable log10 scaling for depth plots"
    )

    add_modkit_args(parser)
    add_motif_args(parser)

    return parser.parse_args()
