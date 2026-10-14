import numpy as np
import pandas as pd
import pyranges as pr
from bedrock.config import BEDMETHYL_COLUMNS, DROP_BEDMETHYL_COLUMNS

def load_chr_map(chr_list_file):
    df = pd.read_csv(chr_list_file, sep="\t", header=None, names=["ref_chr", "assigned_chr"], lineterminator="\n")
    df["ref_chr"] = df["ref_chr"].astype(str).str.strip()
    # one-column list (samtools idxstats | awk '{print $1}'): keep the reference names
    df["assigned_chr"] = df["assigned_chr"].fillna(df["ref_chr"]).astype(str).str.strip()
    return dict(zip(df["ref_chr"], df["assigned_chr"]))

def load_fai(fai_file, chr_map):
    df = pd.read_csv(
        fai_file, sep="\t", header=None,
        names=["Chromosome", "length", "offset", "linebases", "linewidth"]
    )
    df["Chromosome"] = df["Chromosome"].replace(chr_map)
    return df

# columns kept from the bedMethyl (the rest are never used) and their compact types:
# text columns become categories (~10x less memory), counts int32, percent float32
BED_KEEP = [c for c in BEDMETHYL_COLUMNS
            if c not in DROP_BEDMETHYL_COLUMNS and c not in ("strand_pos_start", "strand_pos_end")]
BED_CATEGORIES = ("Chromosome", "mod_code", "strand")


def _bedmethyl_sep(path):
    """Tab-only files (modkit >= 0.5) read much faster; older modkit uses spaces after column 10."""
    with open(path) as fh:
        first = fh.readline()
    return "\t" if len(first.rstrip("\n").split("\t")) >= len(BEDMETHYL_COLUMNS) else r"\s+"


def _remap_chromosomes(series, chr_map):
    """Rename chromosomes on a categorical column (cheap: renames categories, not rows)."""
    cats = list(series.cat.categories)
    new = [chr_map.get(str(c), c) for c in cats]
    if len(set(new)) == len(new):
        return series.cat.rename_categories(new)
    return series.astype(str).replace(chr_map).astype("category")   # two names map to one


def read_bedmethyl(path, sample_name, chr_map):
    sep = _bedmethyl_sep(path)
    with open(path) as fh:
        first = fh.readline().rstrip("\n")
    n_found = len(first.split("\t")) if sep == "\t" else len(first.split())
    n_expected = len(BEDMETHYL_COLUMNS)
    if n_found < n_expected:
        raise ValueError(
            f"[ERROR] {sample_name}: BEDMethyl file has {n_found} columns, "
            f"but {n_expected} expected. File format may be incompatible. "
            f"Recheck if in Bedmethyl format (remember only ModKit bedmethyl format is compatible with this BEDROCK version)."
        )
    if n_found > n_expected:
        print(f"[WARN] {sample_name}: BEDMethyl file has {n_found} columns; "
              f"only the first {n_expected} will be used.")

    idx = [BEDMETHYL_COLUMNS.index(c) for c in BED_KEEP]
    dtypes = {c: ("category" if c in BED_CATEGORIES else "float32" if c == "percent_mod" else "int32")
              for c in BED_KEEP}
    bed = pd.read_csv(path, sep=sep, header=None, usecols=idx, quotechar="'", quoting=0,
                      names=BEDMETHYL_COLUMNS + [f"extra_{i}" for i in range(max(0, n_found - n_expected))],
                      dtype=dtypes)
    bed = bed[BED_KEEP]
    bed["Chromosome"] = _remap_chromosomes(bed["Chromosome"], chr_map)
    return bed


def read_depth(path, chr_map):
    depth = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                        names=["Chromosome", "Position", "Coverage_depth"],
                        dtype={"Chromosome": "category", "Position": "int32", "Coverage_depth": "int32"})
    depth["Chromosome"] = _remap_chromosomes(depth["Chromosome"], chr_map)
    return depth


def load_samples(spreadsheet, chr_map):
    meta = pd.read_csv(
        spreadsheet, sep="\t", header=None,
        names=["sample_name", "sample_type", "bedmethyl_path", "depth_path"]
    )
    samples = {}
    print(f"[INFO] bedMethyl columns used: {', '.join(BED_KEEP)}")

    for _, row in meta.iterrows():
        sample_name = str(row["sample_name"]).strip().strip("'\"")
        sample_type = str(row["sample_type"]).strip().strip("'\"")

        bed = read_bedmethyl(row.bedmethyl_path, sample_name, chr_map)
        unmapped = bed["Chromosome"].isna().sum()
        if unmapped > 0:
            print(f"[WARN] {sample_name}: {unmapped} BED entries have no chromosome name")

        depth = read_depth(row.depth_path, chr_map)

        samples[sample_name] = {
            "type": sample_type,
            "bed": bed,
            "depth": depth,
            "mean_depth": depth.Coverage_depth.mean()
        }

        print(f"Loaded {sample_name} | {len(bed):,} bedMethyl rows | "
              f"Mean depth {depth.Coverage_depth.mean():.2f}", flush=True)

    return samples


def concat_samples(samples, key):
    """One table for all samples (key = "bed" or "depth"), keeping text columns as categories.

    For "bed" the per-sample tables are released afterwards, so the data is held once.
    """
    from pandas.api.types import union_categoricals

    frames = [s[key] for s in samples.values()]
    names = list(samples)
    cat_cols = [c for c in frames[0].columns if isinstance(frames[0][c].dtype, pd.CategoricalDtype)]
    for c in cat_cols:
        cats = union_categoricals([f[c] for f in frames], ignore_order=True).categories
        for f in frames:
            f[c] = f[c].cat.set_categories(cats)
    out = pd.concat(
        [f.assign(sample_name=pd.Categorical.from_codes(np.full(len(f), i, dtype="int32"),
                                                        categories=names))
         for i, f in enumerate(frames)],
        ignore_index=True,
    )
    if key == "bed":
        for s in samples.values():
            s["bed"] = None
    return out


def load_genes_as_pyranges(gff_path):
    df = pd.read_csv(
        gff_path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "Chromosome", "source", "feature",
            "start", "end", "score",
            "strand", "frame", "attributes"
        ]
    )

    df = df[df["feature"] == "gene"].copy()

    if df.empty:
        raise ValueError("No 'gene' features found in GFF file")

    df["gene_id"] = df["attributes"].str.extract(r"ID=([^;]+)")

    if df["gene_id"].isna().any():
        n_missing = df["gene_id"].isna().sum()
        print(
            f"[WARN] {n_missing} gene entries missing gene_id; "
            "loaded but without name"
        )

    df["Start"] = df["start"] - 1
    df["End"] = df["end"]

    df = df[[
        "Chromosome",
        "Start",
        "End",
        "strand",
        "gene_id"
    ]]

    return pr.PyRanges(df)