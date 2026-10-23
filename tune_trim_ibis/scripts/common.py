import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# external inputs
MASTER_DB = "/local/storage/all_birds_highcov_lowcov/Finch_Database_V2.csv"
BIG_BIRDS = "/local/storage/BigBirds/BigBirds.csv"
VCF_SAMPLES = "/local/storage/all_birds_highcov_lowcov/samples.txt"

# the founder of the Big Bird lineage
FOUNDER = "04Dap5110"
DESCENDANT_LINEAGES = ["A line", "B line", "Big Birds"]

# column names of the IBIS .seg output
SEG_COLUMNS = ["sample1", "sample2", "chrom", "phys_start_pos", "phys_end_pos", "IBD_type", "genetic_start_pos", "genetic_end_pos",
               "genetic_seg_length", "marker_count", "error_count", "error_density"]

# plot colors and a theme more like R
COLOR_MAP = {"fortis": "#1C86EE", "scandens": "#43CD80", "fuliginosa": "#C6E2FF", "hybrid": "#AB82FF", "A line": "#FF8C00", "founder": "#FF3030",
             "magnirostris": "#EEB422", "conirostris": "#8B4513", "B line": "#FFAEB9", "unknown": "#808080", "Big Birds": "#800080"}
PLOT_THEME = {"font.family": "DejaVu Sans", "axes.facecolor": "white", "axes.edgecolor": "#DCDCDC", "axes.linewidth": 1.2, "axes.grid": True,
              "grid.color": "#EAEAEA", "grid.linestyle": "-", "grid.linewidth": 0.5, "xtick.color": "black", "ytick.color": "black", "xtick.major.size": 0,
              "ytick.major.size": 0}


def set_plot_theme():
    sns.set_style("whitegrid", rc=PLOT_THEME)
    plt.rcParams.update(PLOT_THEME)


def load_master():
    """Daphne birds from the master database with missing first years filled and the newest lineage classifications."""
    master = pd.read_csv(MASTER_DB)
    master.drop_duplicates(subset=["meaningful.unique"], inplace=True)
    master = master[master["Island"] == "Daphne"].copy()

    # fill in missing years with the minimum year, year born, or last year
    master["First.year"] = master["First.year"].fillna(master["First.year.min"])
    master["First.year"] = master["First.year"].fillna(master["Year.born"])
    master["First.year"] = master["First.year"].fillna(master["Last.year"])
    master["Last.year"] = master["Last.year"].fillna(master["First.year"])

    # update the newest lineage classifications and add generations from the founder
    bbs_lines = pd.read_csv(BIG_BIRDS)
    bbs_lines.drop_duplicates(subset=["meaningful.unique"], inplace=True)
    bbs_lines = bbs_lines.set_index("meaningful.unique")
    master["lineage"] = master["meaningful.unique"].map(bbs_lines["lineage"])
    master["lineage"] = master["lineage"].fillna(master["Species"])
    master = master.drop(columns=["Species"]).rename(columns={"lineage": "Species"})
    master["Generation"] = master["meaningful.unique"].map(bbs_lines["Generation"])
    return master


def first_offspring_year(master):
    """First year any descendant of the founder was observed (1987)."""
    return master[master["Species"].isin(DESCENDANT_LINEAGES)]["First.year"].min()


def negative_controls(master):
    """Birds first observed before any of the founder's descendants (excluding the founder)."""
    is_early = master["First.year"] < first_offspring_year(master)
    return master[is_early & (master["meaningful.unique"] != FOUNDER)]


def positive_controls(master):
    """The founder's first generation offspring."""
    return master[master["Generation"] == 1]


def load_founder_ibis(prefix):
    """Read IBIS .seg and .coef files, keeping only pairs with the founder and adding the partner as "ind"."""
    chunks = pd.read_csv(f"{prefix}.seg", sep="\t", header=None, names=SEG_COLUMNS, chunksize=5_000_000)
    seg = pd.concat(chunk[(chunk["sample1"] == FOUNDER) | (chunk["sample2"] == FOUNDER)] for chunk in chunks)
    seg["ind"] = np.where(seg["sample1"] == FOUNDER, seg["sample2"], seg["sample1"])

    coef = pd.read_csv(f"{prefix}.coef", sep="\t")
    coef = coef[(coef["Individual1"] == FOUNDER) | (coef["Individual2"] == FOUNDER)].copy()
    coef["ind"] = np.where(coef["Individual1"] == FOUNDER, coef["Individual2"], coef["Individual1"])
    return seg, coef


def read_bim(prefix):
    """Read a plink .bim file (chromosomes are the integer codes from make_map.py)."""
    return pd.read_csv(f"{prefix}.bim", sep="\t", header=None, names=["chrom", "snp", "cm", "pos", "a0", "a1"])


def count_snps(segs, bim, chrom="chrom", start="phys_start_pos", end="phys_end_pos"):
    """Number of SNPs in the bim file between start and end (inclusive) for each segment."""
    counts = np.zeros(len(segs), dtype=int)
    positions = {c: np.sort(group["pos"].to_numpy()) for c, group in bim.groupby("chrom")}
    starts, ends = segs[start].to_numpy(), segs[end].to_numpy()
    for c, idx in segs.groupby(chrom).indices.items():
        pos = positions.get(c, np.array([]))
        counts[idx] = np.searchsorted(pos, ends[idx], side="right") - np.searchsorted(pos, starts[idx], side="left")
    return counts


def load_chr_map(path="dataset/chr_map.txt"):
    """Map from the integer chromosome codes used by IBIS to chromosome names (without "chr")."""
    chr_map = pd.read_csv(path, sep=" ", header=None, names=["Key", "Value"])
    chr_map["Key"] = chr_map["Key"].str.replace("chr", "")
    return chr_map.set_index("Value")["Key"].to_dict()


def assign_y_levels(chrom_df, n_ind=None):
    """Stack segments on one chromosome into non-overlapping tracks for plotting (optionally scaled by n_ind)."""
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
    if n_ind:
        chrom_df["y"] /= n_ind
    return chrom_df
