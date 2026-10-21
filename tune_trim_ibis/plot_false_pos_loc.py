import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# make plot theme more like R
custom_params = {"font.family": "DejaVu Sans", "axes.facecolor": "white", "axes.edgecolor": "#DCDCDC", "axes.linewidth": 1.2, "axes.grid": True,
    "grid.color": "#EAEAEA", "grid.linestyle": "-", "grid.linewidth": 0.5, "xtick.color": "black", "ytick.color": "black", "xtick.major.size": 0,
    "ytick.major.size": 0}
sns.set_style("whitegrid", rc=custom_params)
plt.rcParams.update(custom_params)

# first get the master database and filter
master = pd.read_csv("/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv")
master.drop_duplicates(subset=["meaningful.unique"], inplace=True)
master = master[master["Island"] == "Daphne"].copy()

# fill in missing years with the minimum year, year born, or last year
master["First.year"] = master["First.year"].fillna(master["First.year.min"])
master["First.year"] = master["First.year"].fillna(master["Year.born"])
master["First.year"] = master["First.year"].fillna(master["Last.year"])
master["Last.year"] = master["Last.year"].fillna(master["First.year"])

# get the data and convert positions to Mb
seg = pd.read_pickle("output/seg_all_founder.pkl")
seg["phys_start_pos"] /= 1e6
seg["phys_end_pos"] /= 1e6
seg["length_Mb"] = seg["phys_end_pos"] - seg["phys_start_pos"]
seg["Year"] = seg["ind"].map(master.set_index("meaningful.unique")["First.year"])

# get the chromosome mapping
chr_map = pd.read_csv("chr_map.txt", sep=" ", header=None, names=["Key", "Value"])
chr_map["Key"] = chr_map["Key"].str.replace("chr", "")
chr_dict = chr_map.set_index("Value")["Key"].to_dict()

# make a function for assigning y-levels to avoid overlaps
def assign_y_levels(chrom_df):
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
    return chrom_df

# setup the plot
fig, axes = plt.subplots(3, 6, figsize=(22, 9), sharex=False, sharey=False)
axes = axes.flatten()
ax = 0

# store high false-positive regions
false_pos_seg = pd.DataFrame(columns=["chrom", "phys_start_pos", "phys_end_pos"])

neg_con_seg = seg[seg["Year"] < 1987].copy()
for chrom in [c for c in chr_dict if c in neg_con_seg["chrom"].unique()]:
    # plot the locations of the false positive IBD segments
    neg_con_seg_chr = assign_y_levels(neg_con_seg[neg_con_seg["chrom"] == chrom].copy())

    # determine high false-positive regions
    neg_con_seg_chr["diff"] = neg_con_seg_chr["y"].diff().fillna(1)
    neg_con_seg_chr["seg_id"] = (neg_con_seg_chr["diff"] != 1).cumsum()
    pivot_table = neg_con_seg_chr.groupby("seg_id").agg({"y": "count", "phys_start_pos": "min", "phys_end_pos": "max"}).reset_index()

    # Groups containing at least 6 segments
    high_fp_ids = pivot_table.loc[pivot_table["y"] >= 6, "seg_id"]
    pivot_table = pivot_table[pivot_table["y"] >= 6]

    # add to false positive segments
    pivot_table["chrom"] = chrom
    pivot_table = pivot_table[["chrom", "phys_start_pos", "phys_end_pos"]]
    if false_pos_seg.empty:
        false_pos_seg = pivot_table
    else:
        false_pos_seg = pd.concat([false_pos_seg, pivot_table], ignore_index=True)

    # Color segments according to their group
    for _, row in neg_con_seg_chr.iterrows():
        color = "red" if row["seg_id"] in high_fp_ids.values else "steelblue"
        axes[ax].hlines(y=row["y"], xmin=row["phys_start_pos"], xmax=row["phys_end_pos"], color=color, lw=1.5)

    axes[ax].set_xlabel(f"Physical Position (Mb) on Chromosome {chr_dict[chrom]}")
    axes[ax].set_ylabel("IBD Segment Count")
    sns.despine(top=True, right=True, left=True, bottom=True, ax=axes[ax])
    ax += 1

axes[17].set_visible(False)
plt.tight_layout()
plt.savefig("trimming_plots/false_pos_loc.png", dpi=450, bbox_inches="tight", pad_inches=0.05)
plt.close()

# create a bed file for the false positive segments
false_pos_seg.sort_values(by=["chrom", "phys_start_pos", "phys_end_pos"], inplace=True)
false_pos_seg["chrom"] = false_pos_seg["chrom"].apply(lambda x: f"chr{x}")
false_pos_seg["phys_start_pos"] = (false_pos_seg["phys_start_pos"] * 1e6).astype(int)
false_pos_seg["phys_end_pos"] = (false_pos_seg["phys_end_pos"] * 1e6).astype(int)
false_pos_seg.to_csv("trimming/exclude.bed", sep="\t", index=False, header=False)

# create a bed file for the original segments
seg = pd.read_pickle("output/seg_founder_neg_max_error.pkl")
seg = seg[["chrom", "phys_start_pos", "phys_end_pos"]].copy()
seg.sort_values(by=["chrom", "phys_start_pos", "phys_end_pos"], inplace=True)
seg["chrom"] = seg["chrom"].apply(lambda x: f"chr{x}")
seg["name"] = seg.index
seg.to_csv("trimming/features.bed", sep="\t", index=False, header=False)