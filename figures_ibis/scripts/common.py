import importlib.util
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from pandas.api.types import CategoricalDtype
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# paths are found from this file's location, so scripts work from any working directory
FIGURES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(FIGURES)
IBIS = f"{REPO}/tune_trim_ibis"
DATA = f"{FIGURES}/data"

# reuse the bird metadata and plot theme from the IBIS pipeline (loaded by path since both files are named common.py)
_spec = importlib.util.spec_from_file_location("ibis_common", f"{IBIS}/scripts/common.py")
sys.modules["ibis_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["ibis_common"])
from ibis_common import COLOR_MAP, load_chr_map, load_master, read_bim, set_plot_theme

# axis label for the fraction of SNPs in IBD with the founder
IBD_LABEL = "Genome-wide Proportion of SNPs in IBD"

# groups of birds shown in the figures
FOCAL_SPECIES = ["fortis", "scandens"]
GROUPS = ["fortis", "scandens", "A line", "B line", "Big Birds"]

# true species epithets get italicized in figure legends; lineage groups (A line, B line,
# Big Birds, founder, hybrid, unknown) are not species and stay upright
SPECIES_NAMES = {"fortis", "scandens", "conirostris", "magnirostris", "fuliginosa"}


############################ IBIS data ############################
def load_bim():
    return read_bim(f"{IBIS}/plink/dataset")


def load_chromosomes():
    """Map from the integer chromosome codes used by IBIS to chromosome names (without "chr")."""
    return load_chr_map(f"{IBIS}/dataset/chr_map.txt")


def load_trimmed_segments():
    """Trimmed IBD segments shared with the founder (positions in Mb)."""
    return pd.read_pickle(f"{IBIS}/output/trim_seg_founder_neg_max_error.pkl")


def load_founder_coef():
    """Kinship coefficients between the founder and every other bird in the IBIS run."""
    return pd.read_pickle(f"{IBIS}/output/coef_all_founder.pkl")


def load_all_coef():
    """Kinship coefficients between all pairs of birds in the IBIS run."""
    return pd.read_csv(f"{IBIS}/output/finches.coef", sep="\t")


def add_bird_info(df, master, col="ind"):
    """Add species and first/last year from the master database, dropping birds without years."""
    info = master.set_index("meaningful.unique")
    df["Species"] = df[col].map(info["Species"])
    df["Year"] = df[col].map(info["First.year"])
    df["Last_Year"] = df[col].map(info["Last.year"])
    return df.dropna(subset=["Year", "Last_Year"])


def expand_years_alive(df):
    """One row per year each bird was alive (Year_alive)."""
    return df.assign(Year_alive=df.apply(lambda row: range(int(row.Year), int(row.Last_Year) + 1), axis=1)).explode("Year_alive")


def founder_ibd_per_bird(master):
    """Every bird in the IBIS run with its fraction of SNPs in (trimmed) IBD with the founder (0 if none)."""
    snp_counts = load_trimmed_segments().groupby("ind")["SNP_count"].sum()
    birds = load_founder_coef()[["ind"]].copy()
    birds["fraction_snps_ibd"] = (birds["ind"].map(snp_counts) / len(load_bim())).fillna(0)
    return add_bird_info(birds, master)


def yearly_summary(birds):
    """Per year and species: mean fraction of SNPs in IBD, birds alive, birds with IBD, and relative abundance."""
    alive = expand_years_alive(birds)
    grouped = alive.groupby(["Year_alive", "Species"])["fraction_snps_ibd"]
    summary = pd.DataFrame({"fraction_snps_ibd_mean": grouped.mean(), "count_all": grouped.count(),
                            "count_ibd": grouped.apply(lambda x: (x > 0).sum())}).reset_index().rename(columns={"Year_alive": "Year"})
    summary["ibd_fraction"] = summary["count_ibd"] / summary["count_all"]
    summary["spp_fraction"] = summary["count_all"] / summary.groupby("Year")["count_all"].transform("sum")
    return summary


def order_groups(df, groups, by="Species"):
    """Keep only the given species/lineages, in that order (for plotting and legends)."""
    df = df[df["Species"].isin(groups)].copy()
    df["Species"] = df["Species"].astype(CategoricalDtype(categories=groups, ordered=True))
    return df.sort_values(by=by)


def add_bill_pcs(df):
    """PCs of standardized bill length, depth, and width, signed so larger values mean bigger / more pointed beaks."""
    bill = df[["bill.length", "bill.depth", "bill.width"]]
    pcs = PCA(n_components=2).fit_transform(StandardScaler().fit_transform(bill))
    df["PC(Beak Size)"] = pcs[:, 0] * np.sign(np.corrcoef(pcs[:, 0], bill.sum(axis=1))[0, 1])
    df["PC(Beak Pointedness)"] = pcs[:, 1] * np.sign(np.corrcoef(pcs[:, 1], bill["bill.length"] / bill["bill.depth"])[0, 1])
    return df


############################ IBD along the genome ############################
def chromosome_lengths(bim):
    """Position of the last SNP on each chromosome (Mb)."""
    return bim.groupby("chrom")["pos"].max() / 1e6


def compute_freq_profile(df, species, chrom, chrom_len_mb, n_ind, bin_size=0.1):
    """Fraction of birds of a species with an IBD segment overlapping each bin along a chromosome."""
    bins = np.arange(0, chrom_len_mb + bin_size, bin_size)
    bin_centers = bins[:-1] + bin_size / 2
    counts = np.zeros(len(bins) - 1)
    sub = df[(df["Species"] == species) & (df["chrom"] == chrom)]
    if sub.empty:
        return bin_centers, counts
    for ind, grp in sub.groupby("ind"):
        covered = np.zeros(len(bins) - 1, dtype=bool)
        for _, row in grp.iterrows():
            lo = max((np.searchsorted(bins, row["phys_start_pos"], side="right") - 1), 0)
            hi = min((np.searchsorted(bins, row["phys_end_pos"], side="left")), len(covered))
            covered[lo:hi] = True
        counts += covered.astype(float)
    freq = counts / n_ind[species]
    return bin_centers, freq


def ibd_frequency_profiles(master, year, chrom_lengths, bin_size=0.1):
    """IBD frequency profiles for fortis and scandens alive in a given year, keyed by (species, chromosome)."""
    def alive_in(df):
        df = expand_years_alive(add_bird_info(df, master))
        return df[df["Species"].isin(FOCAL_SPECIES) & (df["Year_alive"] == year)].copy()

    segs = alive_in(load_trimmed_segments())
    n_ind = alive_in(load_founder_coef()[["ind"]].copy())["Species"].value_counts()
    return {(species, chrom): compute_freq_profile(segs, species, chrom, chrom_lengths[chrom], n_ind, bin_size)
            for species in FOCAL_SPECIES for chrom in chrom_lengths.index}


def set_genome_plot_theme():
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.facecolor": "white", "axes.edgecolor": "none", "axes.grid": False,
                         "xtick.color": "black", "ytick.color": "black", "xtick.major.size": 0, "ytick.major.size": 0})


def add_freq_colorbar(fig, rect, cmap, vmax, label="IBD segment frequency"):
    """Vertical colorbar of IBD segment frequency from 0 to vmax."""
    sm = cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=0, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=fig.add_axes(rect))
    cb.ax.tick_params(labelsize=13)
    tick_vals = np.linspace(0, vmax, 5)
    cb.set_ticks(tick_vals)
    cb.set_ticklabels([f"{v:.2f}" for v in tick_vals])
    cb.set_label(label, fontsize=14, labelpad=10)
    return cb


############################ plot formatting ############################
def italicize_species_legend(leg):
    for text in leg.get_texts():
        if text.get_text() in SPECIES_NAMES:
            text.set_fontstyle("italic")


def style_legend(ax, **kwargs):
    """White legend with a thin frame and italicized species names."""
    leg = ax.legend(facecolor="white", framealpha=1.0, **kwargs)
    leg.get_frame().set_linewidth(0.3)
    italicize_species_legend(leg)
    return leg


def save_figure(path, fig=None):
    """Save (and close) a figure with the settings shared by all figures (path relative to figures_ibis/, e.g. main/x.png)."""
    fig = fig or plt.gcf()
    fig.savefig(f"{FIGURES}/{path}", dpi=450, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def degree_cutoffs(coef, degrees):
    """Kinship coefficients halfway between each degree of relatedness and the next (closest first)."""
    return [(coef[coef["Degree"] == deg]["Kinship_Coefficient"].min() + coef[coef["Degree"] == deg + 1]["Kinship_Coefficient"].max()) / 2 for deg in degrees]


def format_year_axis(ax):
    ax.set_xlabel(None)
    ax.set_xticks(np.arange(1975, 2025, 5))
    ax.set_xlim(1978, 2024)


# load the drought survival list (Erik's published table) for Daphne fortis, with IBD and beak size
def load_drought_fortis(master):
    drought_surv = pd.read_csv(f"{DATA}/table_of_drought_fortis.csv").set_index("meaningful.unique")["drought_surv"]
    birds = master[master["meaningful.unique"].isin(drought_surv.index)].copy()
    birds["outcome"] = birds["meaningful.unique"].map(drought_surv)
    birds = birds[birds["Species"] == "fortis"].copy()
    ibd = founder_ibd_per_bird(master).set_index("ind")["fraction_snps_ibd"]
    birds["fraction_snps_ibd"] = birds["meaningful.unique"].map(ibd).fillna(0)
    birds["has_ibd"] = birds["fraction_snps_ibd"] > 0
    birds["beak_sum"] = birds["bill.length"] + birds["bill.depth"] + birds["bill.width"]
    return birds
