import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.cluster.hierarchy import leaves_list, linkage
from common import GROUPS, PEAK_REGIONS, WORK, donor_alleles, load_birds, load_genotypes, save_figure, set_plot_theme

# heatmap of every bird's genotype (copies of the donor allele) at the SNPs in the FST peaks, the rows ordered by a tree of all the birds (UPGMA:
# average linkage, cityblock distance on the plotted SNPs), drawn at the side with each bird's group
SNPS = "diagnostic"  # "diagnostic": SNPs telling the donors (conirostris + magnirostris) and scandens (no IBD) apart, or "pruned": the LD-pruned SNPs
ROW_GROUPS = ["magnirostris_dap", "conirostris_esp", "founder", "scandens_ibd", "scandens_no_ibd"]
FOUNDER_ROWS = 2  # the founder's single row drawn this many rows tall, so it can be seen
GAP = 2  # blank columns between peaks
SHOWN_SCANDENS = 150  # scandens without IBD shown (random, fixed seed), so the few conirostris and magnirostris rows stay readable; all are used for
SEED = 1              # choosing the diagnostic SNPs
COLORS = ["#ffffcc", "#fd8d3c", "#bd0026"]  # 0, 1, 2 copies

set_plot_theme()
pos, ids, genotypes = load_genotypes()
group = load_birds().set_index("ind")["group"].loc[ids].to_numpy()
is_donor, diagnostic = donor_alleles(genotypes, group)
mag = np.where(is_donor, genotypes, 2 - genotypes)
scandens = np.flatnonzero(group == "scandens_no_ibd")
hidden = np.setdiff1d(scandens, np.random.default_rng(SEED).choice(scandens, min(SHOWN_SCANDENS, len(scandens)), replace=False))
group = group.astype(object)
group[hidden] = "hidden"  # not plotted

# the chosen SNPs in each peak, in position order
peaks = pd.read_csv(PEAK_REGIONS, sep="\t")
if SNPS == "diagnostic":
    keep = diagnostic
else:
    keep = np.isin(pos, pd.read_csv(f"{WORK}/pruned.prune.in", header=None)[0].str.split(":").str[1].astype(int))
peak_cols = [np.flatnonzero(keep & (pos >= p["start"]) & (pos < p["end"])) for _, p in peaks.iterrows()]
print(f"{SNPS} SNPs per peak:", [len(c) for c in peak_cols])


def with_gaps(matrix):
    """The peaks' columns side by side, separated by GAP blank columns."""
    blocks = []
    for k, cols in enumerate(peak_cols):
        if k:
            blocks.append(np.full((matrix.shape[0], GAP), np.nan))
        blocks.append(matrix[:, cols])
    return np.hstack(blocks)


def fill_missing(rows):
    """Missing genotypes filled with the column mean, for clustering."""
    return np.where(np.isnan(rows), np.nanmean(rows, axis=0), rows)


def draw_tree(ax, tree, leaf_y):
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


# rows ordered by the tree (the founder's row drawn FOUNDER_ROWS rows tall)
shown = np.flatnonzero(np.isin(group, ROW_GROUPS))
tree = linkage(fill_missing(mag[shown][:, np.concatenate(peak_cols)]), method="average", metric="cityblock")
order = shown[leaves_list(tree)]
rows = np.concatenate([np.repeat(r, FOUNDER_ROWS if group[r] == "founder" else 1) for r in order])

fig, (ax_t, ax_g, ax) = plt.subplots(1, 3, figsize=(15.5, 13), gridspec_kw={"width_ratios": [0.12, 0.025, 1], "wspace": 0.01})
ax.imshow(np.ma.masked_invalid(with_gaps(mag[rows])), aspect="auto", interpolation="nearest", cmap=ListedColormap(COLORS), vmin=-0.5, vmax=2.5)

# the tree, its leaves at the middle of each bird's rows
leaf_y = {r: np.flatnonzero(rows == r).mean() for r in order}
draw_tree(ax_t, tree, [leaf_y[r] for r in shown])
ax_t.set_ylim(len(rows) - 0.5, -0.5)

# group color bar (the founder outlined in black)
for k, r in enumerate(rows):
    ax_g.axhspan(k - 0.5, k + 0.5, color=GROUPS.loc[group[r], "color"], linewidth=0)
k = np.flatnonzero(group[rows] == "founder")[0]
ax_g.add_patch(plt.Rectangle((0.04, k - 0.5), 0.92, FOUNDER_ROWS, fill=False, edgecolor="black", linewidth=1, zorder=3))
ax_g.set_ylim(len(rows) - 0.5, -0.5)
ax_g.set_xlim(0, 1)
ax_g.axis("off")

# peak labels at the top
start = 0
for (_, p), cols in zip(peaks.iterrows(), peak_cols):
    ax.text(start + (len(cols) - 1) / 2, -1, f"{p['start'] / 1e6:.2f}-{p['end'] / 1e6:.2f} Mb", ha="center", va="bottom", fontsize=12)
    start += len(cols) + GAP
ax.set_xticks([])
ax.set_yticks([])
ax.grid(False)
for spine in ax.spines.values():
    spine.set_visible(False)

# legends: genotypes, and groups
ax.legend([plt.Rectangle((0, 0), 1, 1, color=c) for c in COLORS], ["0", "1", "2"], title="Copies of the\nconirostris /\nmagnirostris allele",
          loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, fontsize=11, title_fontsize=11)
fig.legend([plt.Rectangle((0, 0), 1, 1, facecolor=GROUPS.loc[g, "color"], edgecolor="black" if g == "founder" else "none") for g in ROW_GROUPS],
           [GROUPS.loc[g, "label"] for g in ROW_GROUPS], loc="upper left", bbox_to_anchor=(1.01, 0.75), bbox_transform=ax.transAxes, frameon=False,
           fontsize=11)
ax.set_title(f"Genotypes in the chr. 6 FST peak{'s' if len(peaks) > 1 else ''} ({sum(len(c) for c in peak_cols)} {SNPS} SNPs; {SHOWN_SCANDENS} of the "
             f"{len(scandens)} scandens without IBD shown)", pad=28)
save_figure("peak_heatmap.png", fig)
