import pandas as pd
from common import load_master, negative_controls, read_bim

# convert IBIS segments (1-based, inclusive) to a sorted bed file (0-based, half-open)
def to_bed(segs, name, path):
    bed = pd.DataFrame({"chrom": "chr" + segs["chrom"].astype(str), "start": segs["phys_start_pos"] - 1, "end": segs["phys_end_pos"], "name": name})
    bed.sort_values(by=["chrom", "start", "end"]).to_csv(path, sep="\t", index=False, header=False)

# get all segments shared with the founder and the negative controls in this IBIS run
seg = pd.read_pickle("output/seg_all_founder.pkl")
coef = pd.read_pickle("output/coef_all_founder.pkl")
neg_ids = sorted(set(negative_controls(load_master())["meaningful.unique"]) & set(coef["ind"]))
with open("trimming/neg_control_ids.txt", "w") as f:
    f.writelines(f"{ind}\n" for ind in neg_ids)

# check that no bird has overlapping segments, so coverage counts individuals rather than segments
neg_seg = seg[seg["ind"].isin(neg_ids)].sort_values(by=["ind", "chrom", "phys_start_pos"])
prev_end = neg_seg.groupby(["ind", "chrom"])["phys_end_pos"].shift()
assert not (neg_seg["phys_start_pos"] <= prev_end).any()
to_bed(neg_seg, neg_seg["ind"], "trimming/neg_control.bed")

# the segments to trim (birds sharing more IBD with the founder than any negative control)
seg = pd.read_pickle("output/seg_founder_neg_max_error.pkl")
to_bed(seg, seg.index, "trimming/features.bed")

# chromosome lengths for bedtools genomecov (up to the last SNP)
bim = read_bim("plink/dataset")
genome = bim.groupby("chrom")["pos"].max().reset_index()
genome["chrom"] = "chr" + genome["chrom"].astype(str)
genome.sort_values(by="chrom").to_csv("trimming/genome.txt", sep="\t", index=False, header=False)
