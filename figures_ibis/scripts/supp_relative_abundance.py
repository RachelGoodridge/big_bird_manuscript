import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common import (COLOR_MAP, GROUPS, format_year_axis, founder_ibd_per_bird, load_master,
                    order_groups, save_figure, set_plot_theme, style_legend, yearly_summary)

set_plot_theme()

# relative abundance of each species / lineage per year
master = load_master()
summary = order_groups(yearly_summary(founder_ibd_per_bird(master)), GROUPS, by=["Year", "Species"])

# add Big Birds end point in 2022
bb_extinct = pd.DataFrame({"Year": [2022], "Species": pd.Categorical(["Big Birds"], categories=GROUPS, ordered=True), "spp_fraction": [0.0]})
summary = pd.concat([summary, bb_extinct], ignore_index=True).sort_values(by=["Year", "Species"])


############################
fig, ax = plt.subplots(1, 1, figsize=(8.2, 12.8))
sns.lineplot(x="Year", y="spp_fraction", hue="Species", palette=COLOR_MAP, data=summary, linewidth=1.5, ax=ax, zorder=2)
sns.scatterplot(x="Year", y="spp_fraction", hue="Species", palette=COLOR_MAP, edgecolor="none", data=summary, s=25, ax=ax, zorder=3, legend=False)

ax.set_ylabel("Relative Species Abundance")
format_year_axis(ax)
ax.set_xlim(1978, 2024)
ax.set_ylim(-0.02, 0.88)

# mark the founder's arrival (1981) as a star sitting on the x-axis
ax.plot(1981, 0, marker="*", markersize=12, color=COLOR_MAP["founder"], clip_on=False, zorder=4)

style_legend(ax, loc="upper left")
sns.despine(top=True, right=True, left=False, bottom=False, ax=ax)

# inset: zoom on A line, B line, and Big Birds ([left, bottom, width, height] in axes fraction)
ax_inset = ax.inset_axes([0.65, 0.85, 0.35, 0.35])
for s in ["A line", "B line", "Big Birds"]:
    d = summary[summary["Species"] == s].sort_values("Year")
    ax_inset.plot(d["Year"].values, d["spp_fraction"].values, color=COLOR_MAP[s], linewidth=1.2, zorder=2)
    ax_inset.scatter(d["Year"].values, d["spp_fraction"].values, color=COLOR_MAP[s], s=10, zorder=3)

ax_inset.set_xticks(np.arange(1980, 2025, 10))
ax_inset.set_xlim(1986, 2024)
ax_inset.set_yticks(np.arange(0, 0.08, 0.02))
ax_inset.set_ylim(-0.002, 0.065)
ax_inset.tick_params(labelsize=10, length=0, colors="black")
ax_inset.set_facecolor("white")
ax_inset.set_axisbelow(True)
ax_inset.grid(True, color="#EAEAEA", linestyle="-", linewidth=0.5)
for spine in ax_inset.spines.values():
    spine.set_edgecolor("#DCDCDC")
    spine.set_linewidth(1.2)
ax_inset.spines["top"].set_visible(False)
ax_inset.spines["right"].set_visible(False)

plt.tight_layout(rect=[0, 0, 1, 0.68])
save_figure("supplement/year_vs_relative_abundance.png")
