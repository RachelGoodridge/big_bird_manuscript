import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pandas.api.types import CategoricalDtype
from common import COLOR_MAP, FOCAL_SPECIES, degree_cutoffs, load_all_coef, load_master, save_figure, set_plot_theme, style_legend

set_plot_theme()
descendants_list = ["A line", "B line", "Big Birds"]

# kinship between all pairs, with first year and species of both birds
master = load_master().set_index("meaningful.unique")
coef = load_all_coef()
for i in [1, 2]:
    coef[f"year{i}"] = coef[f"Individual{i}"].map(master["First.year"])
    coef[f"spp{i}"] = coef[f"Individual{i}"].map(master["Species"])

# include both orientations of each pair, then keep fortis/scandens paired with older descendants of the founder
swapped = coef.rename(columns={"Individual1": "Individual2", "Individual2": "Individual1", "year1": "year2", "year2": "year1", "spp1": "spp2", "spp2": "spp1"})
pairs = pd.concat([coef, swapped], ignore_index=True)
pairs = pairs[pairs["spp1"].isin(FOCAL_SPECIES) & pairs["spp2"].isin(descendants_list) & (pairs["year1"] > pairs["year2"])]

# closest relative in each lineage for each bird, and which lineage is closest
pivot = pairs.pivot_table(index=["Individual1", "spp1", "year1"], columns="spp2", values="Kinship_Coefficient", aggfunc="max").reset_index()
pivot[descendants_list] = pivot[descendants_list].fillna(0)
pivot["max_line"] = pivot[descendants_list].idxmax(axis=1).astype(CategoricalDtype(categories=descendants_list, ordered=True))
pivot["max_kin"] = pivot[descendants_list].max(axis=1)
pivot = pivot.sort_values(by="max_line")
pivot["year1"] = pivot["year1"].astype(int)

# kinship cutoffs between degrees of relatedness (high kinship → low kinship order; 5th degree and beyond are shown together)
cutoffs = degree_cutoffs(coef, range(1, 5))


############################
# plot max line kinship, colored by which line is closer
fig, axes = plt.subplots(figsize=(12, 5), ncols=2, sharex=True, sharey=True)
for ax, spp in zip(axes, FOCAL_SPECIES):
    sns.scatterplot(x="year1", y="max_kin", hue="max_line", data=pivot[pivot["spp1"] == spp], palette=COLOR_MAP, s=50, alpha=0.6, ax=ax, zorder=2)

    # circle the first degree Big Birds
    target_point = pivot[(pivot["spp1"] == spp) & (pivot["max_line"] == "Big Birds") & (pivot["max_kin"] > cutoffs[0])]
    if not target_point.empty:
        sns.scatterplot(data=target_point, x="year1", y="max_kin", facecolor="none", edgecolor="black", s=50, linewidth=1, zorder=3)

    ax.set_xlabel(None)
    ax.set_ylabel("Kinship coefficient to closest relative")
    ax.set_title(spp, fontstyle="italic")
    style_legend(ax, title=None, loc="upper left")
    ax.grid(axis="y")
    ax.set_xticks(np.arange(1975, 2025, 5))
    ax.set_ylim(-0.002, 0.25)
    ax.set_xlim(1978, 2024)

    # show the cutoffs between degrees of relatedness
    for c in cutoffs:
        ax.axhline(y=c, color="lightgray", linewidth=1, zorder=1)

axes[0].get_legend().remove()

# add degree-of-relatedness labels on right side of second panel
y_min, y_max = -0.002, 0.25
boundaries = [y_max] + cutoffs + [y_min]
midpoints = [(boundaries[i] + boundaries[i + 1]) / 2 for i in range(len(boundaries) - 1)]
ax2 = axes[1].twinx()
ax2.set_ylim(y_min, y_max)
ax2.set_yticks(midpoints)
ax2.set_yticklabels(["1st", "2nd", "3rd", "4th", "5th+"])
ax2.tick_params(axis="y", length=0)
ax2.grid(False)

sns.despine(top=True, right=True)
plt.tight_layout()
save_figure("main/year_vs_max_kin.png")
