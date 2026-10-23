import pandas as pd
from common import FOUNDER, SEG_COLUMNS, count_snps, load_master, negative_controls, read_bed, read_bim, subtract_regions

# the founder and each negative control act as a "focal" bird, analyzed the same way as the founder
coef = pd.read_pickle("output/coef_all_founder.pkl")
neg_ids = set(negative_controls(load_master())["meaningful.unique"]) & set(coef["ind"])
focal_birds = neg_ids | {FOUNDER}

# segments between each focal bird and every other bird (pairs of two focal birds are kept once for each)
cols = ["sample1", "sample2", "chrom", "phys_start_pos", "phys_end_pos"]
pairs = []
for chunk in pd.read_csv("output/finches.seg", sep="\t", header=None, names=SEG_COLUMNS, usecols=cols, chunksize=5_000_000):
    for focal, partner in [("sample1", "sample2"), ("sample2", "sample1")]:
        pairs.append(chunk[chunk[focal].isin(focal_birds)].rename(columns={focal: "focal", partner: "partner"}))
pairs = pd.concat(pairs, ignore_index=True).astype({"focal": "category", "partner": "category"})

# fraction of SNPs each pair shares in IBD
bim = read_bim("plink/dataset")
pairs["SNP_count"] = count_snps(pairs, bim)
fraction_snps_ibd = pairs.groupby(["focal", "partner"], observed=True)["SNP_count"].transform("sum") / len(bim)

# keep pairs sharing more IBD than the max negative control error for the founder (same threshold as plot_ibis.py)
founder_frac = pd.read_pickle("output/seg_all_founder.pkl").groupby("ind")["SNP_count"].sum() / len(bim)
neg_max_error = founder_frac.reindex(sorted(neg_ids)).fillna(0).max()
pairs = pairs[fraction_snps_ibd > neg_max_error].reset_index(drop=True)
print(f"negative control max error: {neg_max_error:.5f}, kept {len(pairs)} segments for {pairs['focal'].nunique()} focal birds")

# trim the same false positive regions as for the founder (see trim_segments.py)
trim_seg = subtract_regions(pairs, read_bed("trimming/exclude.bed"))

# recount SNPs, dropping pieces with 1 or fewer SNPs, and convert positions to Mb (as in trim_seg_founder_neg_max_error.pkl)
trim_seg["SNP_count"] = count_snps(trim_seg, bim)
trim_seg = trim_seg[trim_seg["SNP_count"] > 1].copy()
trim_seg["phys_start_pos"] /= 1e6
trim_seg["phys_end_pos"] /= 1e6
trim_seg.reset_index(drop=True).to_pickle("output/trim_seg_controls.pkl")
