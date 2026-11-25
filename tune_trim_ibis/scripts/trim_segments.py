import pandas as pd
from common import count_snps, database_sex, high_coverage_regions, load_chr_map, load_master, negative_controls, read_bim, subtract_regions, write_bed

# get all segments shared with the founder and the negative controls in this IBIS run
seg = pd.read_pickle("output/seg_all_founder.pkl")
coef = pd.read_pickle("output/coef_all_founder.pkl")
neg_ids = set(negative_controls(load_master())["meaningful.unique"]) & set(coef["ind"])
neg_seg = seg[seg["ind"].isin(neg_ids)]

# check that no bird has overlapping segments, so coverage counts individuals rather than segments
sorted_seg = neg_seg.sort_values(by=["ind", "chrom", "phys_start_pos"])
assert not (sorted_seg["phys_start_pos"] <= sorted_seg.groupby(["ind", "chrom"])["phys_end_pos"].shift()).any()

# exclude regions where more than 5% of negative controls share IBD with the founder (likely false positives). On chrZ only the male negative
# controls count (5% of the males): female chrZ is hemizygous and much of it is missing, so the males are the comparable controls there
chrz = next(code for code, name in load_chr_map().items() if name == "Z")
master = load_master()
male_ids = set(master.loc[database_sex(master) == 0, "meaningful.unique"]) & neg_ids
autosome_exclude = high_coverage_regions(neg_seg[neg_seg["chrom"] != chrz], 0.05 * len(neg_ids))
chrz_exclude = high_coverage_regions(neg_seg[(neg_seg["chrom"] == chrz) & neg_seg["ind"].isin(male_ids)], 0.05 * len(male_ids))
exclude = pd.concat([autosome_exclude, chrz_exclude], ignore_index=True)
write_bed(exclude, "trimming/exclude.bed")
print(f"excluding {len(autosome_exclude)} autosomal regions shared with more than {0.05 * len(neg_ids):g} of {len(neg_ids)} negative controls")
print(f"excluding {len(chrz_exclude)} chrZ regions shared with more than {0.05 * len(male_ids):g} of {len(male_ids)} male negative controls")

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
