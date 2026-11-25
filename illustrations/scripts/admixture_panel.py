import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory, offset_copy
from illustration_common import COLOR_MAP, REPO, tint
from ibis_common import load_master

# genome-wide ADMIXTURE ancestry (K = 9) of each bird as stacked bars, one block per population (figure 1C); the same plot as
# ancestry_estimates/whole_data_new/admixture/02_plot_admixture.R, drawn with matplotlib so that it can share the other panels' text. Drawn only as part of figure 1 (fig_1_combined.py)

# ---------------- settings you may want to change ----------------
ADMIXTURE = f"{REPO}/ancestry_estimates/whole_data_new/admixture/finches"
K = 9
# the ancestry cluster (column of the .Q file, from 0) each named population is dominated by, and its color; every other cluster is gray. The
# conirostris cluster is the founder's red, slightly lightened (as in the R plot)
SPECIES_CLUSTER = {"fortis": 0, "scandens": 5, "conirostris": 1}
CLUSTER_COLOR = {"fortis": COLOR_MAP["fortis"], "scandens": COLOR_MAP["scandens"], "conirostris": tint(COLOR_MAP["founder"], 0.13)}
OTHER_COLOR = "#CCCCCC"
# populations left to right: (species, island); fortis and scandens only birds alive before YEAR_CUTOFF (before any offspring of the founder)
POPULATIONS = [("scandens", "Daphne"), ("fortis", "Daphne"), ("conirostris", "Espanola"), ("founder", "Daphne"), ("A line", "Daphne"),
               ("B line", "Daphne"), ("Big Birds", "Daphne")]
YEAR_CUTOFF = 1987
GAP = 2                       # blank bars between populations
GAP_AFTER = {"conirostris": 1}  # the founder (itself a conirostris) sits almost against the conirostris
ISLAND_NAMES = {"Espanola": "Española"}
ISLAND_LINE_BELOW = 2.0       # the line over each island's name, this many label sizes below the populations' names
ISLAND_SIZE = 0.8             # the islands' names, as a fraction of the label size
ISLAND_NAME_GAP = 6           # points between each island's line and its name
LEGEND_SPACING = dict(columnspacing=0.85, handletextpad=0.4)  # between the legend's entries, and between each box and its text (in font sizes)
# ------------------------------------------------------------------


def load_ancestry():
    """Ancestry proportions (columns 0 .. K - 1) of the birds shown, with their population."""
    q = pd.read_csv(f"{ADMIXTURE}.{K}.Q", sep=r"\s+", header=None)
    q["ind"] = pd.read_csv(f"{ADMIXTURE}.fam", sep=r"\s+", header=None)[1]
    master = load_master(island=None)[["meaningful.unique", "Species", "Island", "year"]]
    data = q.merge(master, left_on="ind", right_on="meaningful.unique")
    shown = pd.DataFrame(POPULATIONS, columns=["Species", "Island"])
    data = data.merge(shown, on=["Species", "Island"])
    return data[~data["Species"].isin(["fortis", "scandens"]) | (data["year"] < YEAR_CUTOFF)]


def draw_admixture(ax, label_size=15, legend=True):
    """The stacked bars on one axis, with the populations' names underneath and a legend above."""
    data = load_ancestry()
    clusters = list(range(K))
    start = 0
    extents = {}  # each population's first and last x
    for species, _ in POPULATIONS:
        birds = data[data["Species"] == species]
        means = birds[clusters].mean()
        birds = birds.sort_values(means.idxmax(), ascending=False)  # birds ordered by the population's main ancestry

        # stacked from the bottom: the named clusters from the population's most to least common, then all the others together (gray)
        named = sorted(SPECIES_CLUSTER, key=lambda s: -means[SPECIES_CLUSTER[s]])
        layers = [(birds[SPECIES_CLUSTER[s]].to_numpy(), CLUSTER_COLOR[s]) for s in named]
        layers.append((1 - sum(values for values, _ in layers), OTHER_COLOR))
        edges = start + np.arange(len(birds) + 1)
        # each layer is drawn from 0 up to its top, the tallest first, so that no thin gaps show between the layers
        tops = np.cumsum([values for values, _ in layers], axis=0)
        for top, (_, color) in list(zip(tops, layers))[::-1]:
            ax.stairs(top, edges, baseline=0, fill=True, color=color, linewidth=0)

        extents[species] = (start, start + len(birds))
        middle = start + len(birds) / 2
        if species == "founder":
            ax.text(middle, 1.0, "*", ha="center", va="bottom", fontsize=label_size, fontweight="bold", color="black")  # the single founder
        else:
            named_species = species in SPECIES_CLUSTER  # species names in italics, with the genus
            ax.text(middle, -0.07, f"G. {species}" if named_species else species, ha="center", va="top", fontsize=label_size, color="black",
                    style="italic" if named_species else "normal")
        start += len(birds) + GAP_AFTER.get(species, GAP)

    # the island of each run of neighboring populations (the founder with the Daphne ones after it): a line under their names, the island under
    # it. Each line spans exactly its populations (neighboring islands' lines are apart by the gap between their populations)
    names_bottom = offset_copy(blended_transform_factory(ax.transData, ax.transAxes), fig=ax.figure, y=-ISLAND_LINE_BELOW * label_size,
                               units="points")
    runs = []
    for species, island in POPULATIONS:
        if runs and runs[-1][0] == island:
            runs[-1][2] = species
        else:
            runs.append([island, species, species])
    for island, first, last in runs:
        x0, x1 = extents[first][0], extents[last][1]
        ax.plot([x0, x1], [-0.07, -0.07], transform=names_bottom, color="black", linewidth=1.2, solid_capstyle="butt", clip_on=False)
        ax.text((x0 + x1) / 2, -0.07, ISLAND_NAMES.get(island, island), transform=offset_copy(names_bottom, fig=ax.figure, y=-ISLAND_NAME_GAP, units="points"),
                ha="center", va="top", fontsize=label_size * ISLAND_SIZE, color="black")

    ax.set_xlim(0, start - GAP_AFTER.get(POPULATIONS[-1][0], GAP))
    ax.set_ylim(0, 1)
    ax.axis("off")
    if legend:
        handles = [Patch(color=CLUSTER_COLOR[s]) for s in SPECIES_CLUSTER] + [Patch(color=OTHER_COLOR)]
        leg = ax.legend(handles, [f"G. {s}" for s in SPECIES_CLUSTER] + ["other ancestries"], loc="lower left", bbox_to_anchor=(0, 1.08), ncol=4,
                        frameon=False, fontsize=label_size * 0.9, handlelength=1, handleheight=1, borderaxespad=0, **LEGEND_SPACING)
        for text in leg.get_texts()[:len(SPECIES_CLUSTER)]:
            text.set_style("italic")
