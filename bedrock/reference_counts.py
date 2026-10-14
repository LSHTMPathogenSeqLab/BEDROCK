# reference base counts (C+G and A+T, both strands)
import pandas as pd
from Bio import SeqIO


def count_reference_bases(fasta_path, chrom_map=None):
    chrom_map = chrom_map or {}
    records = []
    genome_totals = {"A": 0, "C": 0}

    for record in SeqIO.parse(fasta_path, "fasta"):
        chrom = chrom_map.get(record.id, record.id)
        seq = str(record.seq).upper()
        a, t, c, g = (seq.count(b) for b in "ATCG")
        genome_totals["A"] += a + t
        genome_totals["C"] += c + g
        records.append({"Chromosome": chrom, "A": a, "T": t, "C": c, "G": g})

    return genome_totals, pd.DataFrame(records)
