import importlib.util
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# paths are found from this file's location, so scripts work from any working directory
ASSOCIATION = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(ASSOCIATION)
DATA = f"{ASSOCIATION}/data"
OUTPUT = f"{ASSOCIATION}/output"
RESULTS = f"{OUTPUT}/ibd_phenotype_association.csv"

# reuse the IBD, trait, window, and LMM helpers from figures_ibis (loaded by path since all are named common.py)
_spec = importlib.util.spec_from_file_location("figures_common", f"{REPO}/figures_ibis/scripts/common.py")
sys.modules["figures_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["figures_common"])
from figures_common import (COLOR_MAP, FOCAL_SPECIES, WINDOW_COVERAGE, add_beak_body_ratio, add_bill_pcs, add_body_pc, eigen_grm, lmm_lrt,
                            load_bim, load_chromosomes, load_founder_coef, load_loco_grm, load_master, load_trimmed_segments, make_windows,
                            set_plot_theme, window_coverage)

# traits tested, and the covariates of each (beak size / pointedness also control for body size, i.e. allometry)
TRAITS = ["PC(Body Size)", "PC(Beak Size)", "PC(Beak Pointedness)", "Beak:Body Ratio"]
COVARIATES = {"PC(Body Size)": ["sex", "loco_ibd"], "PC(Beak Size)": ["sex", "loco_ibd", "PC(Body Size)"],
              "PC(Beak Pointedness)": ["sex", "loco_ibd", "PC(Body Size)"], "Beak:Body Ratio": ["sex", "loco_ibd"]}
MEASUREMENTS = ["weight", "wing", "tarsus", "bill.length", "bill.depth", "bill.width"]

# a window is only tested if its IBD (test) group has at least this many birds, with at least MIN_PER_SEX of each sex
MIN_TEST = {"fortis": 15, "scandens": 10}
MIN_PER_SEX = 3
FDR = 0.05
MANHATTAN_COLORS = {"fortis": ["#206BB5", "#0A3259"], "scandens": ["#3C9E67", "#194D30"]}  # alternating chromosome colors


def load_birds(spp):
    """Birds of a species in the IBIS run with all six measurements and known sex, with the body / beak size PCs and the beak:body ratio fit
    within the species."""
    master = load_master()
    in_ibis = master["meaningful.unique"].isin(load_founder_coef()["ind"])
    birds = master[(master["Species"] == spp) & in_ibis].dropna(subset=MEASUREMENTS + ["sex"]).rename(columns={"meaningful.unique": "ind"})
    return add_beak_body_ratio(add_bill_pcs(add_body_pc(birds.copy()))).reset_index(drop=True)


def _ibd_snps_per_chromosome():
    """SNPs in (trimmed) IBD with the founder for every bird with IBD (bird x IBIS chromosome code), and the SNPs on each chromosome."""
    snps = load_bim().groupby("chrom").size()
    ibd = load_trimmed_segments().groupby(["ind", "chrom"])["SNP_count"].sum().unstack(fill_value=0).reindex(columns=snps.index, fill_value=0)
    return ibd, snps


def loco_ibd_fractions():
    """Fraction of SNPs in (trimmed) IBD with the founder outside each chromosome, for every bird with IBD (bird x IBIS chromosome code)."""
    ibd, snps = _ibd_snps_per_chromosome()
    return (ibd.sum(axis=1).to_numpy()[:, None] - ibd) / (snps.sum() - snps)


REGION_COVERAGE = 0.5  # region test: birds whose IBD covers at least this fraction of the region, vs birds with none in it
REGION_GROUPS = ["No IBD", f"IBD ≥ {REGION_COVERAGE:.0%}"]

# region of chromosome 6 where beak pointedness is associated with IBD in scandens across the whole stretch (windows mostly p < 0.02,
# falling off sharply on either side), tested as a whole in fig_chr6_traits.py and figure 4: chromosome name, start, end (bp)
CHR6_REGION = ("6", 10_000_000, 30_500_000)


def region_test_birds(spp, chrom, start, end):
    """A region-wide version of the window test (chrom: IBIS chromosome code, start / end in bp). The species' tested birds with at least
    REGION_COVERAGE of the region covered by their (summed) IBD segments or none in it (ibd_group; birds with less are left out, like the
    window test's toss group), in the order of the chromosome's LOCO GRM, with the fraction of SNPs in IBD outside the chromosome (loco_ibd),
    and the eigendecomposition of the GRM of these birds."""
    K, ids = load_loco_grm(spp, chrom)
    birds = load_birds(spp).set_index("ind").loc[ids]
    seg = load_trimmed_segments()
    region = pd.DataFrame({"chrom": [chrom], "start": [start], "end": [end]})
    birds["region_ibd"] = window_coverage(seg[seg["ind"].isin(birds.index)], region).set_index("ind")["coverage"].reindex(birds.index).fillna(0)
    birds["ibd_group"] = np.select([birds["region_ibd"] >= REGION_COVERAGE, birds["region_ibd"] == 0], REGION_GROUPS[::-1], default="")
    birds["loco_ibd"] = loco_ibd_fractions()[chrom].reindex(birds.index).fillna(0)
    keep = (birds["ibd_group"] != "").to_numpy()
    return birds[keep], eigen_grm(K[np.ix_(keep, keep)])


def region_lmm(birds, eig, trait):
    """LMM likelihood ratio test of a trait between the birds with and without IBD in the region (birds, eig from region_test_birds), with
    the same GRM and covariates as the window test. Returns beta (the IBD group's difference) and the p-value."""
    has_ibd = (birds["ibd_group"] == REGION_GROUPS[1]).astype(float)
    X = np.column_stack([np.ones(len(birds)), has_ibd] + [birds[c] for c in COVARIATES[trait]])
    return lmm_lrt(birds[trait].to_numpy(dtype=float), X, eig)


def format_lmm(beta, p):
    p_str = f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"
    return f"Linear mixed model: β = {beta:.2f}, {p_str}"


def plot_region_box(ax, birds, trait, beta, p, spp, stats_fontsize="small", point_size=12, seed=1):
    """Every bird's trait value without and with IBD in the region (from region_test_birds), jittered (fixed seed) behind an outlined box plot,
    with beta / p from region_lmm above it."""
    sns.boxplot(data=birds, x="ibd_group", y=trait, order=REGION_GROUPS, fill=False, color="#333333", linewidth=1.5, width=0.55, showfliers=False,
                zorder=2, ax=ax)
    # points drawn after the boxes so the axis limits include every bird
    rng = np.random.default_rng(seed)
    for i, group in enumerate(REGION_GROUPS):
        y = birds.loc[birds["ibd_group"] == group, trait]
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(y)), y, s=point_size, color=COLOR_MAP[spp], alpha=0.5, linewidth=0, zorder=1)
    ax.set_xticks(range(len(REGION_GROUPS)))
    ax.set_xticklabels([f"{group}\n(n = {(birds['ibd_group'] == group).sum()})" for group in REGION_GROUPS])
    ax.set_xlabel(None)
    ax.set_ylabel(trait)
    ax.grid(axis="x", visible=False)
    ax.text(0.5, 1.02, format_lmm(beta, p), transform=ax.transAxes, ha="center", va="bottom", fontsize=stats_fontsize, color="#333333")
    sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)


def genome_layout():
    """Each chromosome's start on a genome-wide axis (chromosomes end to end in genome order) and its length (last SNP), by chromosome name."""
    chr_names = load_chromosomes()
    chroms = list(chr_names.values())
    lengths = load_bim().groupby("chrom")["pos"].max().rename(index=chr_names)[chroms]
    return pd.DataFrame({"start": np.concatenate([[0], np.cumsum(lengths.to_numpy())[:-1]]), "length": lengths.to_numpy()}, index=chroms)


def load_results(genome):
    """Association results with -log10 p, FDR significance, and each window's midpoint on the genome-wide axis."""
    df = pd.read_csv(RESULTS, dtype={"chrom": str})
    df["cum_pos"] = (df["window_start"] + df["window_end"]) / 2 + df["chrom"].map(genome["start"])
    df["neg_log10_p"] = -np.log10(df["p_value"])
    df["significant"] = df["q_value"] < FDR
    return df


def plot_manhattan(ax, sub, genome, spp, y_col, y_label, show_threshold, combine_from="17", point_size=5, tick_fontsize=9):
    """Manhattan plot of one species x trait (sub: its rows of load_results), colored by alternating chromosome. With show_threshold, a line at the
    largest p-value still significant at the FDR cutoff; otherwise significant windows in red and a line at no effect. Chromosomes from
    combine_from on are too small for their own x tick label and share one."""
    palette = {c: MANHATTAN_COLORS[spp][i % 2] for i, c in enumerate(genome.index)}
    if show_threshold:
        sns.scatterplot(data=sub, x="cum_pos", y=y_col, hue="chrom", palette=palette, linewidth=0, s=point_size, alpha=0.7, legend=False, ax=ax, zorder=1)
        if sub["significant"].any():
            ax.axhline(-np.log10(sub.loc[sub["significant"], "p_value"].max()), color="darkgray", linewidth=1, zorder=2)
        ax.set_ylim(bottom=0)
    else:
        sns.scatterplot(data=sub[~sub["significant"]], x="cum_pos", y=y_col, hue="chrom", palette=palette, linewidth=0, s=point_size, alpha=0.7,
                        legend=False, ax=ax, zorder=1)
        ax.scatter(sub.loc[sub["significant"], "cum_pos"], sub.loc[sub["significant"], y_col], color="red", linewidth=0, s=point_size * 1.6, alpha=0.7, zorder=2)
        ax.axhline(0, color="darkgray", linewidth=1, zorder=3)

    # x ticks at each chromosome's center, with the small chromosomes sharing one label
    chroms = list(genome.index)
    split = chroms.index(combine_from)
    end = genome["start"] + genome["length"]
    ax.set_xticks(list(genome["start"].iloc[:split] + genome["length"].iloc[:split] / 2) + [(genome["start"].iloc[split] + end.iloc[-1]) / 2])
    ax.set_xticklabels(chroms[:split] + [f"{combine_from}-28, LGE22"], fontsize=tick_fontsize)
    ax.set_xlim(0, end.iloc[-1])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(y_label)
    ax.grid(False)  # no grid lines (only the threshold / zero line)
    sns.despine(top=True, right=True, ax=ax)


def save_figure(path, fig=None):
    """Save (and close) a figure with the same settings as figures_ibis (path relative to ibd_phenotype_association/figures/)."""
    fig = fig or plt.gcf()
    fig.savefig(f"{ASSOCIATION}/figures/{path}", dpi=450, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
