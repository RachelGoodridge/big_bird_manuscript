import importlib.util
import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import leaves_list, linkage

# paths are found from this file's location, so scripts work from any working directory
HAPLOTYPES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(HAPLOTYPES)
DATA = f"{HAPLOTYPES}/data"
WORK = f"{HAPLOTYPES}/work"  # plink files (large, not tracked)
OUTPUT = f"{HAPLOTYPES}/output"
FIGURES = f"{HAPLOTYPES}/figures"
CANDIDATES = f"{DATA}/candidates.tsv"  # birds passing the IBD and ancestry rules (select_birds.py)
BIRDS = f"{DATA}/birds.tsv"  # every bird analyzed and its group: the candidates passing the regional check (regional_check.py)
FST_GROUPS = f"{DATA}/fst_groups.tsv"  # the two FST groups (fst_groups.py)
FST_WINDOWS = f"{OUTPUT}/fst_windows.tsv"
PEAK_REGIONS = f"{OUTPUT}/fst_peak_regions.tsv"  # FST peaks (fst_vs_recombination.py)
VCF_SAMPLES = "/local/storage/all_birds_highcov_lowcov/samples.txt"  # birds in the phased VCF
GMAP = "/local/storage/gmap_files/chr6.gmap"  # chromosome 6 genetic map (columns pos, chr, cM)

# genome-wide ADMIXTURE ancestry (K = 9) of every bird in the VCF, and the columns of each species' clusters: fortis (V1), scandens (V6), conirostris
# (V2, the Española conirostris and the founder) and magnirostris (V3 + V8, the two magnirostris clusters, summed)
ADMIXTURE = f"{REPO}/ancestry_estimates/whole_data_new/admixture/finches"
ANCESTRY_CLUSTERS = {"fortis": [0], "scandens": [5], "conirostris": [1], "magnirostris": [2, 7]}
MIN_ANCESTRY = 0.5  # birds with less genome-wide ancestry of their own species are left out

# the chromosome 6 region and the IBD rule of the region test in ibd_phenotype_association (figure 4B), loaded by path since all are named common.py
_spec = importlib.util.spec_from_file_location("association_common", f"{REPO}/ibd_phenotype_association/scripts/common.py")
sys.modules["association_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["association_common"])
from association_common import CHR6_REGION, REGION_COVERAGE
from ibis_common import ibd_label
from figures_common import load_chromosomes, load_founder_coef, load_master, load_trimmed_segments, set_plot_theme, species_label, window_coverage

CORE_BIRDS = 100  # reference birds of each species for the regional check (below)

# groups of birds compared (in plotting order): label, color (the project's species colors), island, species (as in load_master), the species whose
# ancestry is checked against MIN_ANCESTRY, and the most birds kept, those with the most ancestry (None: all). The magnirostris are the 39 with the most
# magnirostris ancestry, so that the donors (below) number 50
GROUPS = pd.DataFrame([
    ("fortis",           ibd_label(f"{species_label('fortis')} (No IBD)"), "#1C86EE", "Daphne",   "fortis",       "fortis",       CORE_BIRDS),
    ("scandens_no_ibd",  ibd_label(f"{species_label('scandens')} (No IBD)"), "#43CD80", "Daphne",   "scandens",     "scandens",     None),
    ("scandens_ibd",     ibd_label(f"{species_label('scandens')} (IBD ≥ {REGION_COVERAGE:.0%})"), "#1B6B3A", "Daphne",   "scandens",     "scandens",     None),
    ("conirostris_esp",  f"{species_label('conirostris')} (Española)", "#FF3030", "Espanola", "conirostris",  "conirostris",  None),
    ("founder",          "founder",                                 "#FF3030", "Daphne",   "founder",      "conirostris",  None),
    ("magnirostris_dap", f"{species_label('magnirostris')} (Daphne)",  "#EEB422", "Daphne",   "magnirostris", "magnirostris", 39),
], columns=["group", "label", "color", "island", "species", "ancestry", "max_birds"]).set_index("group")

# the two groups compared (FST scan and diagnostic SNPs): conirostris, the founder, and magnirostris (the donors) vs scandens without IBD. Alleles are
# counted as copies of the donor allele (the allele more common in the donors); SNPs whose frequency differs by at least DIAGNOSTIC between the two
# groups tell the two kinds of haplotype apart
DONORS = ["conirostris_esp", "founder", "magnirostris_dap"]
DIAGNOSTIC = 0.8

# regional check: scandens without IBD left out when, in the chromosome 6 region, they are closer (mean genotype distance, LD-pruned SNPs) to the
# reference fortis than to the reference scandens; the reference birds are the CORE_BIRDS of each with the most ancestry. The fortis are used only
# as this reference (REFERENCE_ONLY), not analyzed
REGIONAL_CHECK = {"scandens_no_ibd": "fortis"}
REFERENCE_ONLY = ["fortis"]

CHROM, START, END = CHR6_REGION


def load_birds():
    return pd.read_csv(BIRDS, sep="\t")


def region_ibd(birds):
    """Fraction of the chromosome 6 region covered by each bird's (summed, trimmed) IBD segments with the founder (0 if none)."""
    chrom_code = {name: code for code, name in load_chromosomes().items()}[CHROM]
    seg = load_trimmed_segments()
    region = pd.DataFrame({"chrom": [chrom_code], "start": [START], "end": [END]})
    coverage = window_coverage(seg[seg["ind"].isin(birds)], region).set_index("ind")["coverage"]
    return coverage.reindex(birds).fillna(0).to_numpy()


def species_ancestry(species):
    """Genome-wide ancestry from a species' clusters (ADMIXTURE, K = 9, summed over ANCESTRY_CLUSTERS) of every bird in the VCF."""
    q = pd.read_csv(f"{ADMIXTURE}.9.Q", sep=r"\s+", header=None)
    ids = pd.read_csv(f"{ADMIXTURE}.fam", sep=r"\s+", header=None)[1]
    return pd.Series(q[ANCESTRY_CLUSTERS[species]].sum(axis=1).to_numpy(), index=ids)


def load_genotypes():
    """Genotypes of every SNP in the region (work/region plink files): SNP positions, bird IDs, and a birds x SNPs array of copies of one of the
    SNP's alleles (0 / 1 / 2, NaN if missing; which allele does not matter, as donor_alleles finds the donor one)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)  # pandas_plink imports pandera the deprecated way
        from pandas_plink import read_plink
        bim, fam, bed = read_plink(f"{WORK}/region", verbose=False)
    return bim["pos"].to_numpy(), fam["iid"].to_numpy(), bed.compute().T


def donor_alleles(genotypes, groups):
    """For each SNP (genotypes: birds x SNPs of allele counts, groups: each bird's group), whether the counted allele is the donor allele, and whether
    the SNP is diagnostic (frequency differing by >= DIAGNOSTIC between the donors and the scandens without IBD, the two FST groups)."""
    freq_donor = np.nanmean(genotypes[np.isin(groups, DONORS)], axis=0) / 2
    freq_sca = np.nanmean(genotypes[groups == "scandens_no_ibd"], axis=0) / 2
    return freq_donor >= freq_sca, np.abs(freq_donor - freq_sca) >= DIAGNOSTIC


def save_figure(path, fig=None):
    """Save (and close) a figure as a png with the same settings as figures_ibis (path relative to chr6_haplotypes/figures/)."""
    fig = fig or plt.gcf()
    fig.savefig(f"{FIGURES}/{path}", dpi=450, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


############################ heatmap of the genotypes in the FST peaks ############################
# every bird's genotype (copies of the donor allele) at the SNPs in the FST peaks, the rows ordered by a tree of all the birds (UPGMA: average
# linkage, cityblock distance on the plotted SNPs), drawn at the side with each bird's group. Drawn by plot_peak_heatmap.py and, as a panel,
# by figures_ibis/scripts/fig_4v2_combined.py
HEATMAP_SNPS = "diagnostic"  # "diagnostic": SNPs telling the donors (conirostris + magnirostris) and scandens (no IBD) apart, or "pruned": the LD-pruned SNPs
HEATMAP_GROUPS = ["magnirostris_dap", "conirostris_esp", "founder", "scandens_ibd", "scandens_no_ibd"]
FOUNDER_ROWS = 2  # the founder's single row drawn this many rows tall, so it can be seen
PEAK_GAP = 2  # blank columns between peaks
SHOWN_SCANDENS = 150  # scandens without IBD shown (random, fixed seed), so the few conirostris and magnirostris rows stay readable; all are used for
SHOWN_SEED = 1        # choosing the diagnostic SNPs
GENOTYPE_COLORS = ["#ffffcc", "#fd8d3c", "#bd0026"]  # 0, 1, 2 copies


def _draw_tree(ax, tree, leaf_y):
    """Dendrogram of a linkage (tree) drawn sideways, its leaves at the rows' positions leaf_y (root on the left)."""
    y = list(leaf_y)
    height = [0.0] * len(leaf_y)
    for a, b, h, _ in tree:
        a, b = int(a), int(b)
        for child in (a, b):
            ax.plot([height[child], h], [y[child], y[child]], color="black", linewidth=0.5)
        ax.plot([h, h], [y[a], y[b]], color="black", linewidth=0.5)
        y.append((y[a] + y[b]) / 2)
        height.append(h)
    ax.set_xlim(tree[-1, 2] * 1.02, 0)
    ax.axis("off")


def plot_peak_heatmap(ax_tree, ax_groups, ax, fontsize=11, peak_fontsize=12, legend_x=1.01, legend_y=(1, 0.75), peak_labels=True, group_labels=None,
                      colors=GENOTYPE_COLORS):
    """The heatmap on three axes side by side: the tree (ax_tree), each row's group color (ax_groups) and the genotypes (ax), with the two
    legends (genotypes, groups) right of the genotypes, their top left corners at legend_x and legend_y (fractions of ax).
    peak_labels: each peak's position written above it; group_labels: {group: label} replacing the groups' usual labels in the legend; colors: of 0, 1 and 2 copies.
    Returns the number of SNPs drawn, of scandens without IBD shown, and of all scandens without IBD."""
    pos, ids, genotypes = load_genotypes()
    group = load_birds().set_index("ind")["group"].loc[ids].to_numpy()
    is_donor, diagnostic = donor_alleles(genotypes, group)
    donor_copies = np.where(is_donor, genotypes, 2 - genotypes)
    scandens = np.flatnonzero(group == "scandens_no_ibd")
    n_shown = min(SHOWN_SCANDENS, len(scandens))
    hidden = np.setdiff1d(scandens, np.random.default_rng(SHOWN_SEED).choice(scandens, n_shown, replace=False))
    group = group.astype(object)
    group[hidden] = "hidden"  # not plotted

    # the chosen SNPs in each peak, in position order
    peaks = pd.read_csv(PEAK_REGIONS, sep="\t")
    if HEATMAP_SNPS == "diagnostic":
        keep = diagnostic
    else:
        keep = np.isin(pos, pd.read_csv(f"{WORK}/pruned.prune.in", header=None)[0].str.split(":").str[1].astype(int))
    peak_cols = [np.flatnonzero(keep & (pos >= p["start"]) & (pos < p["end"])) for _, p in peaks.iterrows()]
    print(f"{HEATMAP_SNPS} SNPs per peak:", [len(c) for c in peak_cols])

    # rows ordered by the tree (missing genotypes filled with the SNP's mean for clustering; the founder's row drawn FOUNDER_ROWS rows tall)
    shown = np.flatnonzero(np.isin(group, HEATMAP_GROUPS))
    clustered = donor_copies[shown][:, np.concatenate(peak_cols)]
    tree = linkage(np.where(np.isnan(clustered), np.nanmean(clustered, axis=0), clustered), method="average", metric="cityblock")
    order = shown[leaves_list(tree)]
    rows = np.concatenate([np.repeat(r, FOUNDER_ROWS if group[r] == "founder" else 1) for r in order])

    # the peaks' columns side by side, separated by PEAK_GAP blank columns
    blocks = []
    for k, cols in enumerate(peak_cols):
        if k:
            blocks.append(np.full((len(rows), PEAK_GAP), np.nan))
        blocks.append(donor_copies[rows][:, cols])
    ax.imshow(np.ma.masked_invalid(np.hstack(blocks)), aspect="auto", interpolation="nearest", cmap=ListedColormap(colors), vmin=-0.5, vmax=2.5)

    # the tree, its leaves at the middle of each bird's rows
    leaf_y = {r: np.flatnonzero(rows == r).mean() for r in order}
    _draw_tree(ax_tree, tree, [leaf_y[r] for r in shown])
    ax_tree.set_ylim(len(rows) - 0.5, -0.5)

    # group color bar (the founder outlined in black)
    for k, r in enumerate(rows):
        ax_groups.axhspan(k - 0.5, k + 0.5, color=GROUPS.loc[group[r], "color"], linewidth=0)
    k = np.flatnonzero(group[rows] == "founder")[0]
    ax_groups.add_patch(plt.Rectangle((0.04, k - 0.5), 0.92, FOUNDER_ROWS, fill=False, edgecolor="black", linewidth=1, zorder=3))
    ax_groups.set_ylim(len(rows) - 0.5, -0.5)
    ax_groups.set_xlim(0, 1)
    ax_groups.axis("off")

    # peak labels at the top
    start = 0
    for (_, p), cols in zip(peaks.iterrows(), peak_cols):
        if peak_labels:
            ax.text(start + (len(cols) - 1) / 2, -1, f"{p['start'] / 1e6:.2f}-{p['end'] / 1e6:.2f} Mb", ha="center", va="bottom", fontsize=peak_fontsize)
        start += len(cols) + PEAK_GAP
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # legends: genotypes, and groups
    ax.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in colors], ["0", "1", "2"],
              title=f"Copies of the\n{species_label('conirostris')} /\n{species_label('magnirostris')} allele",
              loc="upper left", bbox_to_anchor=(legend_x, legend_y[0]), frameon=False, fontsize=fontsize, title_fontsize=fontsize)
    ax.figure.legend([plt.Rectangle((0, 0), 1, 1, facecolor=GROUPS.loc[g, "color"], edgecolor="black" if g == "founder" else "none") for g in HEATMAP_GROUPS],
                     [(group_labels or {}).get(g, GROUPS.loc[g, "label"]) for g in HEATMAP_GROUPS], loc="upper left", bbox_to_anchor=(legend_x, legend_y[1]), bbox_transform=ax.transAxes,
                     frameon=False, fontsize=fontsize)
    return sum(len(c) for c in peak_cols), n_shown, len(scandens)
