import os
import sys
import pandas as pd
from common import (CONTROLS, FOUNDER, IBIS_PREFIX, MERGE_GAP_CM, MIN_LENGTH, ORIGINAL, SEG_COLUMNS, count_snps, filter_and_trim, focal_pairs,
                    load_focal_birds, original_min_length, original_threshold, read_bed, read_bim)

# the founder and each Española conirostris as focal birds, filtered and trimmed as the original run did for its negative controls
# (tune_trim_ibis/scripts/control_segments.py). With --min-length, segments must also be longer than the longest segment the founder shares
# with any negative control, after merging segments separated by small gaps (common.merged_length; applied after the fraction threshold, before
# trimming), with outputs in min_length_merge1cM/ instead of output/
min_length = "--min-length" in sys.argv[1:]
out_dir = MIN_LENGTH if min_length else f"{CONTROLS}/output"
min_cm = original_min_length() if min_length else None
os.makedirs(out_dir, exist_ok=True)

focal_birds = load_focal_birds()
seg = pd.read_csv(f"{IBIS_PREFIX}.seg", sep="\t", header=None, names=SEG_COLUMNS,
                  usecols=["sample1", "sample2", "chrom", "phys_start_pos", "phys_end_pos", "genetic_start_pos", "genetic_end_pos", "genetic_seg_length"])
pairs = focal_pairs(seg, focal_birds)

# the fraction of SNPs each pair shares in IBD (the SNPs are the original run's), with the original run's threshold and false positive regions
bim = read_bim(f"{CONTROLS}/plink/dataset")
threshold = original_threshold()
trim_seg = filter_and_trim(pairs, bim, threshold, read_bed(f"{ORIGINAL}/trimming/exclude.bed"), min_cm)
trim_seg.to_pickle(f"{out_dir}/trim_seg_espanola.pkl")

# summary: for each focal bird, its largest IBD fraction with any bird before filtering, and the birds left after filtering and trimming
fraction = pd.Series(count_snps(pairs, bim), index=pairs.index).groupby([pairs["focal"], pairs["partner"]]).sum() / len(bim)
after = trim_seg.groupby(["focal", "partner"])["SNP_count"].sum() / len(bim)
summary = pd.DataFrame({"max_fraction_before": fraction.groupby(level="focal").max().reindex(focal_birds).fillna(0),
                        "partners_after": after.groupby(level="focal").size().reindex(focal_birds).fillna(0).astype(int),
                        "max_fraction_after": after.groupby(level="focal").max().reindex(focal_birds).fillna(0)})
print(f"threshold (original negative control max error): {threshold:.5f}" + (f", minimum segment length: {min_cm:.3f} cM" if min_length else "") + (f" (gaps of at most {MERGE_GAP_CM:g} cM merged)" if min_length else ""))
print(summary.to_string(float_format=lambda v: f"{v:.5f}"))

# without the length filter, the founder's trimmed segments with Daphne birds should be exactly the original run's
# (tune_trim_ibis/output/trim_seg_founder_neg_max_error.pkl); its pairs with the Española conirostris are new to this run
if not min_length:
    cols = ["chrom", "phys_start_pos", "phys_end_pos", "SNP_count"]
    original = pd.read_pickle(f"{ORIGINAL}/output/trim_seg_founder_neg_max_error.pkl").rename(columns={"ind": "partner"})
    mine = trim_seg[(trim_seg["focal"] == FOUNDER) & ~trim_seg["partner"].isin(focal_birds)]
    sort = lambda df: df[["partner"] + cols].sort_values(["partner"] + cols).reset_index(drop=True)
    assert sort(mine).equals(sort(original)), "the founder's trimmed segments differ from the original run"
    print(f"the founder's {len(mine)} trimmed segments are identical to the original run's")
