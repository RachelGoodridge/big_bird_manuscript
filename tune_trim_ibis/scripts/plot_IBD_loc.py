import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common import assign_y_levels, count_snps, load_chr_map, read_bim, set_plot_theme

set_plot_theme()
chr_dict = load_chr_map()
bim = read_bim("plink/dataset")

# get the data set before trimming
seg_untrim = pd.read_pickle("output/seg_founder_neg_max_error.pkl")
seg_untrim["length_Mb"] = (seg_untrim["phys_end_pos"] - seg_untrim["phys_start_pos"]) / 1e6
num_ind = seg_untrim.ind.nunique()

# get the dataset after trimming (bed back to 1-based positions) and convert positions to Mb
trim_seg = pd.read_csv("trimming/trimmed.bed", sep="\t", header=None, names=["chrom", "phys_start_pos", "phys_end_pos", "name"])
trim_seg["chrom"] = trim_seg["chrom"].str.replace("chr", "").astype(int)
trim_seg["phys_start_pos"] += 1
trim_seg["length_Mb"] = (trim_seg["phys_end_pos"] - trim_seg["phys_start_pos"]) / 1e6
trim_seg["SNP_count"] = count_snps(trim_seg, bim)
trim_seg["phys_start_pos"] /= 1e6
trim_seg["phys_end_pos"] /= 1e6

# setup the plot
fig, axes = plt.subplots(6, 5, figsize=(20, 12), sharex=False, sharey=True)
axes = axes.flatten()

for ax, chrom in zip(axes, chr_dict):
    # plot the locations of the IBD segments
    seg_chr = assign_y_levels(trim_seg[trim_seg["chrom"] == chrom].copy(), num_ind)
    for _, row in seg_chr.iterrows():
        ax.hlines(y=row["y"], xmin=row["phys_start_pos"], xmax=row["phys_end_pos"], color="steelblue", lw=1.5)

    # set the x limits to the span of SNPs on the chromosome
    chrom_pos = bim[bim["chrom"] == chrom]["pos"] / 1e6
    ax.set_xlim(chrom_pos.min(), chrom_pos.max())
    ax.set_xlabel(f"Physical Position (Mb) on Chromosome {chr_dict[chrom]}")
    ax.set_ylabel("Freq. of IBD Segs")
    sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)

plt.tight_layout()
plt.savefig("trimming_plots/IBD_loc_post_trim.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()


############################
# Calculate SNPs per Mb before and after trimming
seg_untrim["SNP_per_Mb"] = seg_untrim["SNP_count"] / seg_untrim["length_Mb"]
trim_seg["SNP_per_Mb"] = trim_seg["SNP_count"] / trim_seg["length_Mb"]

# Print summary statistics
print("SNPs per Mb before trimming:")
print(seg_untrim["SNP_per_Mb"].describe())

trim_seg = trim_seg[trim_seg["SNP_count"]>1].copy()
print("\nSNPs per Mb after trimming:")
print(trim_seg["SNP_per_Mb"].describe())

# Create two side-by-side panels
fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8))

# Before trimming
sns.histplot(data=seg_untrim, x="SNP_per_Mb", ax=axes[0])
axes[0].set_ylabel("Number of IBD Segments")
axes[0].set_title("SNPs per Mb Before Trimming")
axes[0].set_xlabel("SNPs per Mb")
sns.despine(top=True, right=True, left=True, bottom=True, ax=axes[0])

# After trimming
sns.histplot(data=trim_seg, x="SNP_per_Mb", ax=axes[1])
axes[1].set_ylabel("Number of IBD Segments")
axes[1].set_title("SNPs per Mb After Trimming")
axes[1].set_xlabel("SNPs per Mb")
sns.despine(top=True, right=True, left=True, bottom=True, ax=axes[1])

# Adjust layout and save as a single figure
plt.tight_layout()
plt.savefig("trimming_plots/SNPs_per_Mb.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()


############################
# output the trimmed segments
trim_seg["ind"] = trim_seg["name"].map(seg_untrim["ind"])
trim_seg.drop("name", axis=1, inplace=True)
trim_seg.reset_index(drop=True, inplace=True)
trim_seg.to_pickle("output/trim_seg_founder_neg_max_error.pkl")
