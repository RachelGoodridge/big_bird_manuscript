import matplotlib.pyplot as plt
import seaborn as sns
from common import COLOR_MAP, count_snps, load_founder_ibis, load_master, negative_controls, read_bim, set_plot_theme

set_plot_theme()

# load the IBIS segments shared with the founder
seg, coef = load_founder_ibis("output/finches")

# fraction of all SNPs that each bird shares IBD with the founder
bim = read_bim("plink/dataset")
seg["SNP_count"] = count_snps(seg, bim)
fraction_snps_ibd = seg.groupby("ind")["SNP_count"].sum() / len(bim)

# get the birds used in IBIS (birds with no segments share 0% IBD)
master = load_master()
master = master[master["meaningful.unique"].isin(coef["ind"])].copy()
master["fraction_snps_ibd"] = master["meaningful.unique"].map(fraction_snps_ibd).fillna(0)


############################
# get the max error rate for negative controls
neg_max_error = negative_controls(master)["fraction_snps_ibd"].max()

# plot generation vs. fraction of SNPs in IBD
plt.figure(figsize=(6.6, 5))
sns.scatterplot(x="Generation", y="fraction_snps_ibd", hue="Species", palette=COLOR_MAP, data=master[~master["Generation"].isna()], s=50, alpha=0.7)
plt.xlabel("Generations from the Founder")
plt.ylabel("Fraction of all SNPs in IBD")
plt.title("IBD Segments Detected with IBIS")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.legend(title=None)
plt.tight_layout()
plt.savefig("tuning_plots/gen_vs_frac_ibd.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# plot generation vs. distribution of segment lengths
plt.figure(figsize=(6.6, 5))
seg["Generation"] = seg["ind"].map(master.set_index("meaningful.unique")["Generation"])
seg["Species"] = seg["ind"].map(master.set_index("meaningful.unique")["Species"])
sns.boxplot(x="Generation", y="genetic_seg_length", hue="Species", palette=COLOR_MAP, data=seg[~seg["Generation"].isna()])
plt.xlabel("Generations from the Founder")
plt.ylabel("Distribution of IBD Segment Lengths (cM)")
plt.title("IBD Segments Detected with IBIS")
sns.despine(top=True, right=True, left=True, bottom=True)
plt.legend(title=None, loc="upper right")
plt.tight_layout()
plt.savefig("tuning_plots/gen_vs_seg_len_dist.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# plot generation vs. kinship coefficient
plt.figure(figsize=(6.6, 5))
coef["Generation"] = coef["ind"].map(master.set_index("meaningful.unique")["Generation"])
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
inds_to_keep = fraction_snps_ibd[fraction_snps_ibd > neg_max_error].index
seg = seg[seg["ind"].isin(inds_to_keep)].copy()
seg.reset_index(drop=True, inplace=True)
seg.to_pickle("output/seg_founder_neg_max_error.pkl")
