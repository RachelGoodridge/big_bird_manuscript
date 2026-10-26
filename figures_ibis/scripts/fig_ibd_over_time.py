import matplotlib.pyplot as plt
import seaborn as sns
from common import (COLOR_MAP, FOCAL_SPECIES, GROUPS, IBD_LABEL, MEAN_IBD_LABEL, format_year_axis, founder_ibd_per_bird,
                    italicize_species_legend, load_master, order_groups, save_figure, set_plot_theme, yearly_summary)

set_plot_theme()

# fraction of SNPs in IBD with the founder for every bird, and per year averages / carrier frequencies
master = load_master()
birds = founder_ibd_per_bird(master)
summary = order_groups(yearly_summary(birds), FOCAL_SPECIES, by=["Year", "Species"])
birds = order_groups(birds, GROUPS)


############################
# three panels: (a) individual birds, (b) yearly mean, (c) IBD segment carrier frequency
fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(18, 5.2))

# (a) individual scatter points
sns.scatterplot(x="Year", y="fraction_snps_ibd", hue="Species", palette=COLOR_MAP, data=birds, s=35, alpha=0.75, ax=ax_a, zorder=2)
ax_a.set_ylabel(IBD_LABEL)
ax_a.set_ylim(-0.02, 1.0)

# (b) average lines and (c) carrier frequency
for ax, y_col, y_label, y_lim in [(ax_b, "fraction_snps_ibd_mean", MEAN_IBD_LABEL, (-0.002, 0.077)),
                                  (ax_c, "ibd_fraction", "IBD Segment Carrier Frequency", (-0.02, 0.77))]:
    sns.lineplot(x="Year", y=y_col, hue="Species", palette=COLOR_MAP, data=summary, linewidth=1.5, ax=ax, zorder=2)
    sns.scatterplot(x="Year", y=y_col, hue="Species", palette=COLOR_MAP, edgecolor="none", data=summary, s=25, ax=ax, zorder=3, legend=False)
    ax.set_ylabel(y_label)
    ax.set_ylim(*y_lim)

for ax in [ax_a, ax_b, ax_c]:
    format_year_axis(ax)
    ax.get_legend().remove()
    sns.despine(top=True, right=True, left=False, bottom=False, ax=ax)

# one legend for all three panels in a row above them, so it covers no data (panel a has every group)
handles, labels = ax_a.get_legend_handles_labels()
leg = fig.legend(handles, labels, loc="upper center", ncol=len(labels), bbox_to_anchor=(0.5, 0.95), frameon=False)
italicize_species_legend(leg)

plt.tight_layout(rect=(0, 0, 1, 0.92))
save_figure("main/figure2a-c.png")
