import pandas as pd
from common import count_snps, high_coverage_regions, load_master, negative_controls, read_bim, subtract_regions, write_bed

# get all segments shared with the founder and the negative controls in this IBIS run
seg = pd.read_pickle("output/seg_all_founder.pkl")
coef = pd.read_pickle("output/coef_all_founder.pkl")
neg_ids = set(negative_controls(load_master())["meaningful.unique"]) & set(coef["ind"])
neg_seg = seg[seg["ind"].isin(neg_ids)]

# check that no bird has overlapping segments, so coverage counts individuals rather than segments
sorted_seg = neg_seg.sort_values(by=["ind", "chrom", "phys_start_pos"])
assert not (sorted_seg["phys_start_pos"] <= sorted_seg.groupby(["ind", "chrom"])["phys_end_pos"].shift()).any()

# exclude regions where more than 5% of negative controls share IBD with the founder (likely false positives)
exclude = high_coverage_regions(neg_seg, 0.05 * len(neg_ids))
write_bed(exclude, "trimming/exclude.bed")
print(f"excluding {len(exclude)} regions shared with more than {0.05 * len(neg_ids):g} of {len(neg_ids)} negative controls")

# trim the segments of birds sharing more IBD with the founder than any negative control (sorted by chromosome and position)
seg = pd.read_pickle("output/seg_founder_neg_max_error.pkl")
seg = seg.loc[pd.DataFrame({"chrom": "chr" + seg["chrom"].astype(str), "start": seg["phys_start_pos"], "end": seg["phys_end_pos"]}).sort_values(by=["chrom", "start", "end"]).index]
trim_seg = subtract_regions(seg[["chrom", "phys_start_pos", "phys_end_pos", "ind"]].reset_index(drop=True), exclude)

# recount SNPs, dropping pieces with 1 or fewer SNPs, and convert positions to Mb
trim_seg["length_Mb"] = (trim_seg["phys_end_pos"] - trim_seg["phys_start_pos"]) / 1e6
trim_seg["SNP_count"] = count_snps(trim_seg, read_bim("plink/dataset"))
trim_seg["SNP_per_Mb"] = trim_seg["SNP_count"] / trim_seg["length_Mb"]
trim_seg = trim_seg[trim_seg["SNP_count"] > 1].copy()
trim_seg["phys_start_pos"] /= 1e6
trim_seg["phys_end_pos"] /= 1e6
trim_seg = trim_seg[["chrom", "phys_start_pos", "phys_end_pos", "length_Mb", "SNP_count", "SNP_per_Mb", "ind"]]
trim_seg.reset_index(drop=True).to_pickle("output/trim_seg_founder_neg_max_error.pkl")
