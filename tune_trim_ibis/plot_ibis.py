import os
from os.path import join
import pandas as pd
from pandas_plink import read_plink
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="pandera")

# make plot theme more like R
custom_params = {"font.family": "DejaVu Sans", "axes.facecolor": "white", "axes.edgecolor": "#DCDCDC", "axes.linewidth": 1.2, "axes.grid": True,
    "grid.color": "#EAEAEA", "grid.linestyle": "-", "grid.linewidth": 0.5, "xtick.color": "black", "ytick.color": "black", "xtick.major.size": 0,
    "ytick.major.size": 0}
sns.set_style("whitegrid", rc=custom_params)
plt.rcParams.update(custom_params)
color_map_hex = {"fortis": "#1C86EE", "scandens": "#43CD80", "fuliginosa": "#C6E2FF", "hybrid": "#AB82FF", "A line": "#FF8C00", "founder": "#FF3030",
                 "magnirostris": "#EEB422", "conirostris": "#8B4513", "B line": "#FFAEB9", "unknown": "#808080", "Big Birds": "#800080"}

# load the IBIS files
seg = pd.read_csv("output/finches.seg", sep="\t", header=None)
seg.columns = ["sample1", "sample2", "chrom", "phys_start_pos", "phys_end_pos", "IBD_type", "genetic_start_pos", "genetic_end_pos",
               "genetic_seg_length", "marker_count", "error_count", "error_density"]
coef = pd.read_csv("output/finches.coef", sep="\t")

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

# filter for the birds used in ibis
master = master[master["meaningful.unique"].isin(coef.ind.to_list())].copy()
master["Generation"] = master["meaningful.unique"].map(bbs_lines.set_index("meaningful.unique")["Generation"])
master["fraction_snps_ibd"] = master["meaningful.unique"].map(pivot_seg.set_index("ind")["fraction_snps_ibd"]).fillna(0)


############################
# get the max error rate for negative controls
neg_con = master[master["First.year"] < 1987].copy()
neg_max_error = neg_con["fraction_snps_ibd"].max()

# plot generation vs. fraction of SNPs in IBD
plt.figure(figsize=(6.4, 4.8))
sns.scatterplot(x="Generation", y="fraction_snps_ibd", hue="Species", palette=color_map_hex, data=master[~master["Generation"].isna()], s=50, alpha=0.7)
plt.xlabel("Generations from the Founder")
plt.ylabel("Fraction of all SNPs in IBD")
plt.title("IBD Segments Detected with IBIS")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.legend(title=None)
plt.tight_layout()
plt.savefig("tuning_plots/gen_vs_frac_ibd.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# plot generation vs. distribution of segment lengths
plt.figure(figsize=(6.4, 4.8))
seg["Generation"] = seg["ind"].map(bbs_lines.set_index("meaningful.unique")["Generation"])
seg["Species"] = seg["ind"].map(master.set_index("meaningful.unique")["Species"])
sns.boxplot(x="Generation", y="genetic_seg_length", hue="Species", palette=color_map_hex, data=seg[~seg["Generation"].isna()])
plt.xlabel("Generations from the Founder")
plt.ylabel("Distribution of IBD Segment Lengths (cM)")
plt.title("IBD Segments Detected with IBIS")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.legend(title=None, loc="upper right")
plt.tight_layout()
plt.savefig("tuning_plots/gen_vs_seg_len_dist.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# plot generation vs. kinship coefficient
plt.figure(figsize=(6.4, 4.8))
coef["Generation"] = coef["ind"].map(bbs_lines.set_index("meaningful.unique")["Generation"])
sns.scatterplot(x="Generation", y="Kinship_Coefficient", hue="Degree", data=coef[~coef["Generation"].isna()], s=50, alpha=0.7)
plt.xlabel("Generations from the Founder (Actual)")
plt.ylabel("Kinship Coefficient")
plt.title("Kinship Coefficients from IBIS")
plt.legend(title="Generation\n  (Inferred)", loc="lower right")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.tight_layout()
plt.savefig("tuning_plots/gen_vs_kinship.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

############################
# first output all segments shared with the founder
seg.drop(columns=["Generation", "Species"]).to_pickle("output/seg_all_founder.pkl")
coef.drop(columns=["Generation"]).to_pickle("output/coef_all_founder.pkl")

# then filter for those above the negative control max error rate
inds_to_keep = pivot_seg[pivot_seg["fraction_snps_ibd"] > neg_max_error]["ind"].to_list()
seg = seg[seg["ind"].isin(inds_to_keep)].copy()
seg.reset_index(drop=True, inplace=True)
seg.to_pickle("output/seg_founder_neg_max_error.pkl")