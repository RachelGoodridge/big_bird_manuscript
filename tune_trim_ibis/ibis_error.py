import os
from os.path import join
import pandas as pd
from pandas_plink import read_plink
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="pandera")

# load the IBIS files
seg = pd.read_csv("tuning/finches.seg", sep="\t", header=None)
seg.columns = ["sample1", "sample2", "chrom", "phys_start_pos", "phys_end_pos", "IBD_type", "genetic_start_pos", "genetic_end_pos",
               "genetic_seg_length", "marker_count", "error_count", "error_density"]
coef = pd.read_csv("tuning/finches.coef", sep="\t")

# filter for the founder and identify the focal individual
seg = seg[(seg["sample1"] == "04Dap5110") | (seg["sample2"] == "04Dap5110")]
seg["ind"] = np.where((seg["sample1"] == "04Dap5110"), seg["sample2"], seg["sample1"])
coef = coef[(coef["Individual1"] == "04Dap5110") | (coef["Individual2"] == "04Dap5110")]
coef["ind"] = np.where((coef["Individual1"] == "04Dap5110"), coef["Individual2"], coef["Individual1"])


############################
# find the total number of SNPs in the dataset
(bim, fam, bed) = read_plink(join("plink", "dataset*.bed"), verbose=False)
bim["chrom"] = bim["chrom"].astype(int)
num_SNPs = bim.shape[0]

# manually count the number of SNPs in each segment
def SNP_counter(chr, start, end, bim=bim):
    temp_bim = bim[(bim["chrom"]==chr) & (bim["pos"]>=start) & (bim["pos"]<=end)]
    return temp_bim.shape[0]

# figure out the fraction of all SNPs in IBD
seg["SNP_count"] = seg.apply(lambda row: SNP_counter(row["chrom"], row["phys_start_pos"], row["phys_end_pos"]), axis=1)
pivot_seg = pd.pivot_table(seg, index="ind", values=["SNP_count"], aggfunc=["sum"]).reset_index()
pivot_seg["fraction_snps_ibd"] = pivot_seg["sum"]["SNP_count"] / num_SNPs


############################
# first get the master database and filter
master = pd.read_csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv")
master.drop_duplicates(subset=["meaningful.unique"], inplace=True)
master = master[master["Island"] == "Daphne"].copy()

# fill in missing years with the minimum year, year born, or last year
master["First.year"] = master["First.year"].fillna(master["First.year.min"])
master["First.year"] = master["First.year"].fillna(master["Year.born"])
master["First.year"] = master["First.year"].fillna(master["Last.year"])
master["Last.year"] = master["Last.year"].fillna(master["First.year"])

# update the newest lineage classifications
bbs_lines = pd.read_csv("/local/storage/BigBirds/BigBirds.csv")
bbs_lines.drop_duplicates(subset=["meaningful.unique"], inplace=True)
master["lineage"] = master["meaningful.unique"].map(bbs_lines.set_index("meaningful.unique")["lineage"])
master["lineage"] = master["lineage"].fillna(master["Species"])
master = master.drop(columns=["Species"]).rename(columns={"lineage": "Species"})

# filter for the birds for tuning
master = master[master["meaningful.unique"].isin(coef.ind.to_list())].copy()
master["Generation"] = master["meaningful.unique"].map(bbs_lines.set_index("meaningful.unique")["Generation"])
master["fraction_snps_ibd"] = master["meaningful.unique"].map(pivot_seg.set_index("ind")["fraction_snps_ibd"]).fillna(0)


############################
# summarize sharing of first generation birds
pos_con = master[master["Generation"] == 1].copy()
pos_con["error_rate"] = (1 - pos_con["fraction_snps_ibd"])

# summarize sharing of pre-founder birds
neg_con = master[master["Generation"].isna()].copy()

# output the mean and max error rate for negative controls
pos_mean_error = pos_con["error_rate"].mean()
neg_mean_error = neg_con["fraction_snps_ibd"].mean()
neg_max_error = neg_con["fraction_snps_ibd"].max()

print(f"{pos_mean_error},{neg_mean_error},{neg_max_error}")
