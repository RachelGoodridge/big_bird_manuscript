import importlib.util
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pandas.api.types import CategoricalDtype
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# paths are found from this file's location, so scripts work from any working directory
FIGURES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(FIGURES)
IBIS = f"{REPO}/tune_trim_ibis"
GENOTYPING = f"{REPO}/locus_genotyping"  # large-effect beak loci positions and genotype calls

# reuse the bird metadata and plot theme from the IBIS pipeline (loaded by path since both files are named common.py)
_spec = importlib.util.spec_from_file_location("ibis_common", f"{IBIS}/scripts/common.py")
sys.modules["ibis_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["ibis_common"])
from ibis_common import COLOR_MAP, FOUNDER, load_chr_map, load_master, negative_controls, read_bim, set_plot_theme

# axis label for the fraction of SNPs in IBD with the founder
IBD_LABEL = "Genome-wide Proportion of SNPs in IBD"
MEAN_IBD_LABEL = "Mean Proportion of SNPs in IBD"

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
    """Add species and the first/last year alive from the master database, dropping birds without years."""
    info = master.set_index("meaningful.unique")
    df["Species"] = df[col].map(info["Species"])
    df["Year"] = df[col].map(info["year"])
    df["Last_Year"] = df[col].map(info["last.year"])
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


def plot_ibd_over_time(axes, birds, summary, point_scale=1, line_scale=1):
    """Figure 2 panels: (a) each bird's fraction of SNPs in IBD, (b) yearly mean in fortis / scandens, (c) IBD segment carrier frequency.
    Point sizes and line widths are multiplied by point_scale and line_scale.
    Returns the legend handles and labels of panel a (it has every group), for one legend shared by all panels."""
    ax_a, ax_b, ax_c = axes

    # (a) individual scatter points
    sns.scatterplot(x="Year", y="fraction_snps_ibd", hue="Species", palette=COLOR_MAP, data=birds, s=35 * point_scale, alpha=0.75, ax=ax_a, zorder=2)
    ax_a.set_ylabel(IBD_LABEL)
    ax_a.set_ylim(-0.02, 1.0)

    # (b) average lines and (c) carrier frequency
    for ax, y_col, y_label, y_lim in [(ax_b, "fraction_snps_ibd_mean", MEAN_IBD_LABEL, (-0.002, 0.077)),
                                      (ax_c, "ibd_fraction", "IBD Segment Carrier Frequency", (-0.02, 0.77))]:
        sns.lineplot(x="Year", y=y_col, hue="Species", palette=COLOR_MAP, data=summary, linewidth=1.5 * line_scale, ax=ax, zorder=2)
        sns.scatterplot(x="Year", y=y_col, hue="Species", palette=COLOR_MAP, edgecolor="none", data=summary, s=25 * point_scale, ax=ax, zorder=3, legend=False)
        ax.set_ylabel(y_label)
        ax.set_ylim(*y_lim)

    for ax in axes:
        format_year_axis(ax)
        ax.tick_params(axis="x", rotation=30)
        ax.get_legend().remove()
        sns.despine(top=True, right=True, left=False, bottom=False, ax=ax)
    return ax_a.get_legend_handles_labels()


def add_bill_pcs(df):
    """PCs of standardized bill length, depth, and width, signed so larger values mean bigger / more pointed beaks."""
    bill = df[["bill.length", "bill.depth", "bill.width"]]
    pcs = PCA(n_components=2).fit_transform(StandardScaler().fit_transform(bill))
    df["PC(Beak Size)"] = pcs[:, 0] * np.sign(np.corrcoef(pcs[:, 0], bill.sum(axis=1))[0, 1])
    df["PC(Beak Pointedness)"] = pcs[:, 1] * np.sign(np.corrcoef(pcs[:, 1], bill["bill.length"] / bill["bill.depth"])[0, 1])
    return df


def add_body_pc(df):
    """First PC of standardized weight, wing, and tarsus, signed so larger values mean heavier birds."""
    pc = PCA(n_components=1).fit_transform(StandardScaler().fit_transform(df[["weight", "wing", "tarsus"]]))[:, 0]
    df["PC(Body Size)"] = pc * np.sign(np.corrcoef(pc, df["weight"])[0, 1])
    return df


def add_beak_body_ratio(df):
    """Beak volume (length x depth x width) relative to body mass, log transformed (closer to normal) and standardized across df (unitless)."""
    beak_volume = df["bill.length"] * df["bill.depth"] * df["bill.width"]
    df["Beak:Body Ratio"] = StandardScaler().fit_transform(np.log(beak_volume / df["weight"]).to_frame()).ravel()
    return df


############################ beak traits over time ############################
def load_measured_fortis_scandens(master):
    """Every fortis / scandens on Daphne with beak length, depth, and width (whether or not in the IBIS run), with years alive."""
    birds = master[master["Species"].isin(FOCAL_SPECIES)].dropna(subset=["bill.length", "bill.depth", "bill.width"]).rename(columns={"meaningful.unique": "ind"})
    return add_bird_info(birds, master)


def yearly_trait_summary(birds, traits, min_birds=15):
    """Per year and species, counting each bird in every year it was alive: mean, sd, coefficient of variation ((sd / mean) * 100,
    as in Grant & Grant) and a point size for the number of birds. Years with fewer than min_birds birds of a species are dropped."""
    alive = expand_years_alive(birds)
    alive = alive[alive.groupby(["Year_alive", "Species"])["ind"].transform("size") >= min_birds]
    grouped = alive.groupby(["Year_alive", "Species"])[list(traits)]
    mean, sd = grouped.mean(), grouped.std()
    n_birds = grouped.size()
    return mean, sd, sd / mean * 100, 5 + 40 * n_birds / n_birds.max()


def plot_trait_means(ax, mean, sd, point_size, trait, point_scale=1, line_scale=1):
    """Yearly mean +/- sd of a trait for fortis and scandens (point size = number of birds).
    Point sizes and line / error bar widths are multiplied by point_scale and line_scale."""
    for spp in FOCAL_SPECIES:
        years = mean.xs(spp, level="Species").index
        m, s = mean.xs(spp, level="Species")[trait], sd.xs(spp, level="Species")[trait]
        ax.errorbar(years, m, yerr=s, color=COLOR_MAP[spp], ecolor=COLOR_MAP[spp], elinewidth=0.8 * line_scale, capsize=2 * line_scale,
                    capthick=plt.rcParams["lines.markeredgewidth"] * line_scale, alpha=0.6, linewidth=0, zorder=1)
        ax.plot(years, m, color=COLOR_MAP[spp], linewidth=1.5 * line_scale, zorder=2, label=spp)
        ax.scatter(years, m, s=point_size.xs(spp, level="Species") * point_scale, color=COLOR_MAP[spp], edgecolor="none", zorder=3)


############################ allometry ############################
def load_allometry_birds():
    """Daphne birds in the IBIS run, the founder, and Española conirostris, with IBD fraction and body / beak size PCs."""
    master = load_master(island=None)
    ibis_birds = set(load_founder_coef()["ind"]) | {FOUNDER}
    daphne = (master["Island"] == "Daphne") & master["meaningful.unique"].isin(ibis_birds)
    conirostris = (master["Island"] == "Espanola") & (master["Species"] == "conirostris")
    birds = master[daphne | conirostris].dropna(subset=["bill.length", "bill.depth", "bill.width", "weight", "wing", "tarsus"]).copy()

    # IBD with the founder (0 for conirostris and the founder, which were not tested)
    ibd = founder_ibd_per_bird(load_master()).set_index("ind")["fraction_snps_ibd"]
    birds["fraction_snps_ibd"] = birds["meaningful.unique"].map(ibd).fillna(0)

    # PCs fit on all of these birds, so both allometry figures share the same axes
    return add_bill_pcs(add_body_pc(birds))


def plot_allometry(birds, fit_birds, fit_label, path):
    """Beak size vs. body size (outline width = fraction of SNPs in IBD) with a regression fit to fit_birds."""
    birds = birds.sort_values(by="fraction_snps_ibd")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(birds["PC(Body Size)"], birds["PC(Beak Size)"], c=birds["Species"].map(COLOR_MAP), alpha=0.5, edgecolors="black", linewidths=birds["fraction_snps_ibd"])

    fit = stats.linregress(fit_birds["PC(Body Size)"], fit_birds["PC(Beak Size)"])
    x = np.sort(fit_birds["PC(Body Size)"])
    ax.plot(x, fit.intercept + fit.slope * x, color="black", linewidth=2)
    ax.text(0.05, 0.99, fit_label, transform=ax.transAxes, ha="left", va="top", fontsize=11.5)
    ax.text(0.05, 0.94, f"$R^2$ = {fit.rvalue ** 2:.3f}", transform=ax.transAxes, ha="left", va="top", fontsize=11.5)

    handles = [Line2D([], [], color=color, marker="o", linestyle="None", markersize=6, label=s) for s, color in COLOR_MAP.items() if s in set(birds["Species"])]
    style_legend(ax, handles=handles, bbox_to_anchor=(1, 1), loc="upper left")
    ax.set_xlabel("PC(Body Size)")
    ax.set_ylabel("PC(Beak Size)")
    sns.despine(top=True, right=True, ax=ax)
    plt.tight_layout()
    save_figure(path, fig)


############################ negative controls as founders ############################
def control_yearly_summary(master):
    """Carrier frequency and mean fraction of SNPs in IBD among fortis / scandens alive each year, for the founder and each
    negative control (each treated as a founder, see tune_trim_ibis/05_control_birds.sh)."""
    # fortis and scandens alive each year in the IBIS run
    birds = add_bird_info(pd.DataFrame({"ind": list(load_founder_coef()["ind"]) + [FOUNDER]}), master)
    alive = expand_years_alive(birds)
    alive = alive[alive["Species"].isin(FOCAL_SPECIES)].drop(columns="Year").rename(columns={"Year_alive": "Year"})
    n_alive = alive.groupby(["Year", "Species"]).size()

    # fraction of SNPs each bird shares in IBD with each focal bird, for every year the bird was alive
    segs = pd.read_pickle(f"{IBIS}/output/trim_seg_controls.pkl")
    pairs = (segs.groupby(["focal", "partner"])["SNP_count"].sum() / len(load_bim())).rename("fraction_snps_ibd").reset_index()
    pairs = expand_years_alive(add_bird_info(pairs, master, col="partner"))
    pairs = pairs[pairs["Species"].isin(FOCAL_SPECIES)].drop(columns="Year").rename(columns={"Year_alive": "Year"})
    sums = pairs.groupby(["focal", "Year", "Species"])["fraction_snps_ibd"].agg(n_ibd="count", ibd_sum="sum")

    # every focal bird x year x species, excluding the focal bird itself from the number alive
    focal_birds = sorted(set(negative_controls(master)["meaningful.unique"]) & set(birds["ind"])) + [FOUNDER]
    summary = pd.MultiIndex.from_product([focal_birds, n_alive.index.get_level_values("Year").unique().sort_values(), FOCAL_SPECIES],
                                         names=["focal", "Year", "Species"]).to_frame(index=False)
    summary = summary.join(sums, on=["focal", "Year", "Species"]).fillna({"n_ibd": 0, "ibd_sum": 0})
    is_self = summary.set_index(["focal", "Year", "Species"]).index.isin(alive.set_index(["ind", "Year", "Species"]).index)
    summary["n_alive"] = summary.join(n_alive.rename("n"), on=["Year", "Species"])["n"].fillna(0) - is_self
    summary = summary[summary["n_alive"] > 0].copy()
    summary["carrier_freq"] = summary["n_ibd"] / summary["n_alive"]
    summary["mean_ibd"] = summary["ibd_sum"] / summary["n_alive"]
    summary["is_founder"] = summary["focal"] == FOUNDER
    return summary


def plot_controls(ax, summary, y_col, y_label):
    """Each negative control (thin), their mean (thick), and the founder (points) over time for fortis and scandens; returns legend handles."""
    colors = {"fortis": {"founder": "#1C86EE", "mean": "#67abed", "control": "#B2CFEC"},
              "scandens": {"founder": "#43CD80", "mean": "#84ddac", "control": "#C6EDD8"}}
    for spp in FOCAL_SPECIES:
        sub = summary[summary["Species"] == spp]
        controls = sub[~sub["is_founder"]]
        for _, bird in controls.groupby("focal"):
            ax.plot(bird["Year"], bird[y_col], color=colors[spp]["control"], linewidth=0.8, alpha=0.8, zorder=1)
        mean = controls.groupby("Year")[y_col].mean()
        ax.plot(mean.index, mean.values, color=colors[spp]["mean"], linewidth=2, zorder=4)
        founder = sub[sub["is_founder"]]
        ax.plot(founder["Year"], founder[y_col], color=colors[spp]["founder"], linewidth=2, marker="o", markersize=3.5, zorder=2)

    handles = [Line2D([0], [0], color=colors[spp][kind], lw=2 if kind != "control" else 0.8, marker="o" if kind == "founder" else None, markersize=3,
                      label=rf"$\it{{{spp}}}$ – {label}")
               for kind, label in [("founder", "founder"), ("mean", "mean of controls"), ("control", "each control")] for spp in FOCAL_SPECIES]
    ax.set_ylabel(y_label)
    format_year_axis(ax)
    ax.set_xlim(1973, 2024)
    ax.set_xticks(range(1975, 2025, 5))
    sns.despine(top=True, right=True, ax=ax)
    return handles


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


# genomic windows for IBD frequencies: non-overlapping, and a bird is in IBD in a window if its segments together cover at least 80% of it
WINDOW_SIZE = 100_000
WINDOW_COVERAGE = 0.8


def make_windows(bim):
    """Non-overlapping windows from the first to the last SNP of each chromosome (the last window may be shorter), with names chrom_index."""
    windows = []
    for chrom, pos in bim.groupby("chrom")["pos"]:
        starts = np.arange(pos.min(), pos.max(), WINDOW_SIZE)
        windows.append(pd.DataFrame({"chrom": chrom, "start": starts, "end": np.minimum(starts + WINDOW_SIZE, pos.max())}))
    windows = pd.concat(windows, ignore_index=True)
    windows["name"] = windows["chrom"].astype(str) + "_" + windows.groupby("chrom").cumcount().astype(str)
    return windows


def ibd_window_frequencies(master):
    """Proportion of fortis / scandens alive each year in IBD with the founder in each window (year x window table per species)."""
    windows = make_windows(load_bim())
    win_starts, win_ends = windows["start"].to_numpy(), windows["end"].to_numpy()

    # trimmed segments of fortis and scandens (positions in bp), and the birds of each species alive each year
    seg = add_bird_info(load_trimmed_segments(), master)
    seg = seg[seg["Species"].isin(FOCAL_SPECIES)].reset_index(drop=True)
    starts_bp, ends_bp = np.rint(seg["phys_start_pos"] * 1e6), np.rint(seg["phys_end_pos"] * 1e6)
    alive = expand_years_alive(founder_ibd_per_bird(master))
    n_alive = alive[alive["Species"].isin(FOCAL_SPECIES)].groupby(["Species", "Year_alive"]).size()

    # length of each segment's overlap with every window it touches
    pieces = []
    for chrom, win in windows.groupby("chrom"):
        chrom_starts, chrom_ends, offset = win["start"].to_numpy(), win["end"].to_numpy(), win.index[0]
        for k in np.flatnonzero(seg["chrom"] == chrom):
            lo, hi = np.searchsorted(chrom_ends, starts_bp[k], side="right"), np.searchsorted(chrom_starts, ends_bp[k], side="left")
            overlap = np.minimum(ends_bp[k], chrom_ends[lo:hi]) - np.maximum(starts_bp[k], chrom_starts[lo:hi])
            pieces.append(pd.DataFrame({"ind": seg.at[k, "ind"], "window": offset + np.arange(lo, hi), "overlap": overlap}))

    # a bird is in IBD in a window if its segments (which never overlap each other) together cover at least 80% of the window
    covered = pd.concat(pieces).groupby(["ind", "window"], as_index=False)["overlap"].sum()
    covered = covered[covered["overlap"] >= WINDOW_COVERAGE * (win_ends - win_starts)[covered["window"]]]

    # count the birds in IBD in each window for every year they were alive
    birds = seg.drop_duplicates("ind").set_index("ind")[["Species", "Year", "Last_Year"]]
    covered_years = expand_years_alive(covered.join(birds, on="ind"))
    years = sorted(covered_years["Year_alive"].unique())
    freqs = {}
    for spp in FOCAL_SPECIES:
        sub = covered_years[covered_years["Species"] == spp]
        counts = np.zeros((len(years), len(windows)))
        np.add.at(counts, (np.searchsorted(years, sub["Year_alive"].astype(int)), sub["window"].to_numpy(dtype=int)), 1)
        freqs[spp] = pd.DataFrame(counts / n_alive[spp].reindex(years).to_numpy()[:, None], index=years, columns=windows["name"])
    return freqs, windows


def set_genome_plot_theme():
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.facecolor": "white", "axes.edgecolor": "none", "axes.grid": False,
                         "xtick.color": "black", "ytick.color": "black", "xtick.major.size": 0, "ytick.major.size": 0})


def add_freq_colorbar(fig, rect, cmap, vmax, label="IBD segment frequency", fontsize=17):
    """Vertical colorbar of IBD segment frequency from 0 to vmax."""
    sm = cm.ScalarMappable(cmap=cmap, norm=Normalize(vmin=0, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=fig.add_axes(rect))
    cb.ax.tick_params(labelsize=fontsize)
    tick_vals = np.linspace(0, vmax, 5)
    cb.set_ticks(tick_vals)
    cb.set_ticklabels([f"{v:.2f}" for v in tick_vals])
    cb.set_label(label, fontsize=fontsize + 2, labelpad=10)
    return cb


def plot_ibd_track(fig, ax, profiles, years, chr_dict, chrom_lengths, vmax, bin_size=0.1, fontsize=15, combine_from="21"):
    """Manhattan-style tracks of IBD segment frequency: chromosomes laid end to end, one track per species and year
    (fortis then scandens, years in the given order), colored by frequency on a 0 to vmax scale (gray where no IBD).
    Chromosomes from combine_from onwards share one x tick label."""
    chroms_ordered = list(chr_dict.keys())
    ax.grid(False)
    ax.set_facecolor("white")

    # gap scales with the largest chromosome so tiny chromosomes stay visually distinct without wasting space on huge ones.
    # Wide enough that the white axes background reads as a clear boundary instead of a thin separator line.
    chrom_gap_mb = chrom_lengths.max() * 0.02
    chrom_offset = {}
    cum_pos = 0.0
    for chrom in chroms_ordered:
        chrom_offset[chrom] = cum_pos
        cum_pos += chrom_lengths[chrom] + chrom_gap_mb
    total_genome_mb = cum_pos - chrom_gap_mb

    track_height = 0.5
    row_pitch = track_height * 1.2  # vertical gap between rows of 0.1 track heights
    track_rows = [(species, yr) for species in FOCAL_SPECIES for yr in years]
    n_rows = len(track_rows)
    cmap = plt.get_cmap("Reds")
    norm = Normalize(vmin=0, vmax=vmax)

    row_y0 = {}
    for row_idx, (species, yr) in enumerate(track_rows):
        y0 = (n_rows - 1 - row_idx) * row_pitch
        row_y0[row_idx] = y0
        for chrom in chroms_ordered:
            x_off = chrom_offset[chrom]
            ax.add_patch(Rectangle((x_off, y0), chrom_lengths[chrom], track_height, facecolor="#d8d8d8", edgecolor="none", zorder=1))
            bc, fq = profiles[yr][(species, chrom)]
            for bcenter, f in zip(bc, fq):
                if f <= 0:
                    continue
                x_lo = x_off + max(bcenter - bin_size / 2, 0.0)
                x_hi = x_off + min(bcenter + bin_size / 2, chrom_lengths[chrom])
                if x_hi - x_lo <= 0:
                    continue
                ax.add_patch(Rectangle((x_lo, y0), x_hi - x_lo, track_height, facecolor=cmap(norm(f)), edgecolor="none", zorder=2))

    # thin outline box around each species' rows (chromosome and row boundaries are shown by the white gaps)
    species_row_indices = {}
    for row_idx, (species, yr) in enumerate(track_rows):
        species_row_indices.setdefault(species, []).append(row_idx)
    for species, rows in species_row_indices.items():
        block_y0 = row_y0[max(rows)]
        block_y1 = row_y0[min(rows)] + track_height
        ax.add_patch(Rectangle((0, block_y0), total_genome_mb, block_y1 - block_y0, facecolor="none", edgecolor="#333333", linewidth=1.0, zorder=5, clip_on=False))

    ax.set_xlim(0, total_genome_mb)
    ax.set_ylim(0, (n_rows - 1) * row_pitch + track_height)

    # the tail of small chromosomes is too narrow for separate labels, so they are combined into one tick
    split_idx = [chr_dict[c] for c in chroms_ordered].index(combine_from)
    individual_chroms = chroms_ordered[:split_idx]
    combined_chroms = chroms_ordered[split_idx:]
    combined_x0 = chrom_offset[combined_chroms[0]]
    combined_x1 = chrom_offset[combined_chroms[-1]] + chrom_lengths[combined_chroms[-1]]
    ax.set_xticks([chrom_offset[c] + chrom_lengths[c] / 2 for c in individual_chroms] + [(combined_x0 + combined_x1) / 2])
    ax.set_xticklabels([str(chr_dict[c]) for c in individual_chroms] + [f"{combine_from}-28, LGE22"], fontsize=fontsize)
    ax.set_xlabel("Chromosome", fontsize=fontsize + 4, labelpad=8)

    # unlabeled minor ticks: one centered in each inter-chromosome gap, plus one at the very start and end of the genome
    boundary_ticks = [0.0] + [(chrom_offset[a] + chrom_lengths[a] + chrom_offset[b]) / 2 for a, b in zip(chroms_ordered[:-1], chroms_ordered[1:])] + [total_genome_mb]
    ax.set_xticks(boundary_ticks, minor=True)
    ax.tick_params(axis="x", which="minor", length=3, color="#aaaaaa", width=0.8)
    ax.tick_params(axis="x", which="major", length=0)
    ax.set_yticks([])

    # two label columns hugging the left edge of the axis, both rotated vertical: species name spanning its rows,
    # and the year of each row closer to the axis (offsets in points so the gap stays fixed regardless of the axes' width)
    row_y_centers = [(n_rows - 1 - i) * row_pitch + track_height / 2 for i in range(n_rows)]
    base_trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
    species_trans = mtransforms.offset_copy(base_trans, fig=fig, x=-30, y=0, units="points")
    year_trans = mtransforms.offset_copy(base_trans, fig=fig, x=-10, y=0, units="points")
    for species, rows in species_row_indices.items():
        y_center = (row_y_centers[min(rows)] + row_y_centers[max(rows)]) / 2
        ax.text(0.0, y_center, r"$\it{" + species + "}$", transform=species_trans, ha="center", va="center", rotation=90, fontsize=fontsize, clip_on=False)
    for row_idx, (species, yr) in enumerate(track_rows):
        ax.text(0.0, row_y_centers[row_idx], f"{yr}", transform=year_trans, ha="center", va="center", rotation=90, fontsize=fontsize, clip_on=False)

    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#aaaaaa")
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(axis="both", length=0)


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
    ax.set_xlim(1977, 2024)
