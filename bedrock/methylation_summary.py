# genome-wide methylation summary (5mC / 5hmC / 6mA)
import pandas as pd
from Bio import SeqIO

MOD_NAME_MAP = {"m": "5mC", "h": "5hmC", "a": "6mA"}
BASE_MOD_CODES = {"C": ["m", "h"], "A": ["a"]}

MOD_COL = "mod"
COV_COL = "mod_score"
POS_COLS = ["Chromosome", "Start_chrom_pos", "strand"]


def methylation_summary(bed_all, outdir="."):
    """Total C / A calls = summed coverage; each modification = summed modified
    calls / total calls for that base."""
    rows = []
    for sample_name, bed in bed_all.groupby("sample_name", observed=True, sort=False):

        for ref_base, codes in BASE_MOD_CODES.items():
            sub = bed[bed["mod_code"].isin(codes)]
            if sub.empty:
                continue

            total_calls = sub.groupby(POS_COLS, observed=True)[COV_COL].max().sum()

            n_by_code = {}
            for code in codes:
                if code not in sub["mod_code"].values:
                    continue
                n_by_code[MOD_NAME_MAP[code]] = sub.loc[sub["mod_code"] == code, MOD_COL].sum()

            if len(n_by_code) > 1:   # C: 5mC + 5hmC combined
                n_by_code["5mC+5hmC combined"] = sum(n_by_code.values())

            for modification, n_mod in n_by_code.items():
                rows.append({
                    "sample_name": sample_name,
                    "ref_base": ref_base,
                    "total_calls": total_calls,
                    "modification": modification,
                    "n_modified": n_mod,
                    "percent_modified": n_mod / total_calls * 100 if total_calls else float("nan"),
                })

    out = pd.DataFrame(rows)
    out.to_csv(f"{outdir}/methylation_summary.tsv", sep="\t", index=False)

    with open(f"{outdir}/methylation_summary.txt", "w") as f:
        f.write("Genome-wide methylation summary (all calls, no mod_score filter)\n")
        f.write("Total calls = summed coverage over all positions called for that base\n")
        f.write("=" * 65 + "\n")
        for sample_name in out["sample_name"].unique():
            f.write(f"\nSample: {sample_name}\n")
            s = out[out["sample_name"] == sample_name]
            for ref_base in s["ref_base"].unique():
                b = s[s["ref_base"] == ref_base]
                f.write(f"  Total {ref_base} calls: {b['total_calls'].iloc[0]:,.0f}\n")
                for _, r in b.iterrows():
                    f.write(
                        f"    {r['modification']:<20s} {r['n_modified']:>12,.0f} modified "
                        f"= {r['percent_modified']:.4f}%\n"
                    )

    print(out)
    return out