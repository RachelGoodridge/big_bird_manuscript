import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common import assign_y_levels, load_chr_map, load_master, negative_controls, read_bed, set_plot_theme

set_plot_theme()
chr_dict = load_chr_map()

# get the negative control segments and the excluded regions (positions in Mb)
seg = pd.read_pickle("output/seg_all_founder.pkl")
neg_seg = seg[seg["ind"].isin(negative_controls(load_master())["meaningful.unique"])][["chrom", "phys_start_pos", "phys_end_pos", "ind"]].copy()
exclude = read_bed("trimming/exclude.bed")
for df in [neg_seg, exclude]:
    df["phys_start_pos"] /= 1e6
    df["phys_end_pos"] /= 1e6
neg_seg["length_Mb"] = neg_seg["phys_end_pos"] - neg_seg["phys_start_pos"]

# setup the plot with one panel per chromosome that has false positive segments
chroms = [c for c in chr_dict if c in neg_seg["chrom"].unique()]
fig, axes = plt.subplots(math.ceil(len(chroms) / 6), 6, figsize=(25.5, 10), sharex=False, sharey=False)
axes = axes.flatten()

for ax, chrom in zip(axes, chroms):
    # plot the locations of the false positive IBD segments
    neg_seg_chr = assign_y_levels(neg_seg[neg_seg["chrom"] == chrom].copy())
    for _, row in neg_seg_chr.iterrows():
        ax.hlines(y=row["y"], xmin=row["phys_start_pos"], xmax=row["phys_end_pos"], color="steelblue", lw=1.5)

    # shade the excluded regions
    for _, row in exclude[exclude["chrom"] == chrom].iterrows():
        ax.axvspan(row["phys_start_pos"], row["phys_end_pos"], color="red", alpha=0.25, lw=0)

    ax.set_xlabel(f"Position (Mb) on Chromosome {chr_dict[chrom]}")
    ax.set_ylabel("IBD Segment Count")
    sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)

for ax in axes[len(chroms):]:
    ax.set_visible(False)

plt.tight_layout()
plt.savefig("trimming_plots/false_pos_loc.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()
