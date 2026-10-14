BEDMETHYL_COLUMNS = [
    "Chromosome",
    "Start_chrom_pos",
    "End_chrom_pos",
    "mod_code",
    "mod_score",
    "strand",
    "strand_pos_start",
    "strand_pos_end",
    "color",
    "mod_cov",
    "percent_mod",
    "mod",
    "canonical_mod",
    "other_mod",
    "delete_mod",
    "fail_mod",
    "variant_mod",
    "no_call",
]

DEPTH_COLUMNS = [
    "Chromosome",
    "Position",
    "Coverage_depth",
]

DROP_BEDMETHYL_COLUMNS = {
    "color",
    "mod_cov",
}


SAMPLE_PALETTE = "colorblind"        # samples
CHROM_PALETTE = "tab20"         # chromosomes
BASE_CALL_PALETTE = "Set2"      # mod / canonical / etc
MODIFICATION_PALETTE = "Dark2"  # 5mC / 5hmC / 6mA
EXTENDED_PALETTE = "hls" #backup
