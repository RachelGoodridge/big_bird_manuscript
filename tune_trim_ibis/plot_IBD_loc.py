import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pandas_plink import read_plink
from os.path import join
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="pandera")

# make plot theme more like R
custom_params = {"font.family": "DejaVu Sans", "axes.facecolor": "white", "axes.edgecolor": "#DCDCDC", "axes.linewidth": 1.2, "axes.grid": True,
    "grid.color": "#EAEAEA", "grid.linestyle": "-", "grid.linewidth": 0.5, "xtick.color": "black", "ytick.color": "black", "xtick.major.size": 0,
    "ytick.major.size": 0}
sns.set_style("whitegrid", rc=custom_params)
plt.rcParams.update(custom_params)

# get the chromosome mapping
chr_map = pd.read_csv("chr_map.txt", sep=" ", header=None, names=["Key", "Value"])
chr_map["Key"] = chr_map["Key"].str.replace("chr", "")
chr_dict = chr_map.set_index("Value")["Key"].to_dict()

# read in the plink files
(bim, fam, bed) = read_plink(join("plink", "dataset*.bed"), verbose=False)
bim["chrom"] = bim["chrom"].astype(int)

# manually count the number of SNPs in each segment
def SNP_counter(chr, start, end, bim):
    temp_bim = bim[(bim["chrom"]==chr) & (bim["pos"]>=start) & (bim["pos"]<=end)]
    return temp_bim.shape[0]

# get the data set before trimming
seg_untrim = pd.read_pickle("output/seg_founder_neg_max_error.pkl")
seg_untrim["length_Mb"] = (seg_untrim["phys_end_pos"] - seg_untrim["phys_start_pos"]) / 1e6
num_ind = seg_untrim.ind.nunique()

# get the dataset after trimming and convert positions to Mb
trim_seg = pd.read_csv("trimming/trimmed.bed", sep="\t", header=None, names=["chrom", "phys_start_pos", "phys_end_pos", "name"])
trim_seg["chrom"] = trim_seg["chrom"].apply(lambda x: int(x.replace("chr", "")))
trim_seg["length_Mb"] = (trim_seg["phys_end_pos"] - trim_seg["phys_start_pos"]) / 1e6
trim_seg["SNP_count"] = trim_seg.apply(lambda row: SNP_counter(row["chrom"], row["phys_start_pos"], row["phys_end_pos"], bim), axis=1)
trim_seg["phys_start_pos"] /= 1e6
trim_seg["phys_end_pos"] /= 1e6

# make a function for assigning y-levels to avoid overlaps
def assign_y_levels(chrom_df, n_ind):
    chrom_df = chrom_df.sort_values(by=["phys_start_pos", "length_Mb"], ascending=[True, False]).reset_index(drop=True)
    tracks = []  # stores the end positions of the last segment in each track
    for idx, row in chrom_df.iterrows():
        placed = False
        for level, last_end in enumerate(tracks):
            if row["phys_start_pos"] > last_end:
                # place segment in this existing track
                tracks[level] = row["phys_end_pos"]
                chrom_df.loc[idx, "y"] = level + 1
                placed = True
                break
        if not placed:
            # need a new track
            tracks.append(row["phys_end_pos"])
            chrom_df.loc[idx, "y"] = len(tracks)
    chrom_df["y"] /= n_ind  # normalize by number of individuals
    return chrom_df

# setup the plot
fig, axes = plt.subplots(6, 5, figsize=(20, 12), sharex=False, sharey=True)
axes = axes.flatten()
ax = 0

for chrom in [key for key in chr_dict.keys()]:
    # plot the locations of the IBD segments
    seg_chr = assign_y_levels(trim_seg[trim_seg["chrom"] == chrom].copy(), num_ind)

    # make the plot
    for _, row in seg_chr.iterrows():
        axes[ax].hlines(y=row["y"], xmin=row["phys_start_pos"], xmax=row["phys_end_pos"], color="steelblue", lw=1.5)

    # figure out the x limits
    xmin = bim[bim["chrom"] == chrom]["pos"].min() / 1e6
    xmax = bim[bim["chrom"] == chrom]["pos"].max() / 1e6

    axes[ax].set_xlim(xmin, xmax)
    axes[ax].set_xlabel(f"Physical Position (Mb) on Chromosome {chr_dict[chrom]}")
    axes[ax].set_ylabel("Freq. of IBD Segs")
    sns.despine(top=True, right=True, left=True, bottom=True, ax=axes[ax])
    ax += 1

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