import matplotlib.pyplot as plt
from common import (FOCAL_SPECIES, GROUPS, founder_ibd_per_bird, italicize_species_legend, load_master, order_groups, plot_ibd_over_time,
                    save_figure, set_plot_theme, yearly_summary)

set_plot_theme()

# fraction of SNPs in IBD with the founder for every bird, and per year averages / carrier frequencies
master = load_master()
birds = founder_ibd_per_bird(master)
summary = order_groups(yearly_summary(birds), FOCAL_SPECIES, by=["Year", "Species"])
birds = order_groups(birds, GROUPS)


############################
# three panels: (a) individual birds, (b) yearly mean, (c) IBD segment carrier frequency
fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))
handles, labels = plot_ibd_over_time(axes, birds, summary)

# one legend for all three panels in a row above them, so it covers no data (panel a has every group)
leg = fig.legend(handles, labels, loc="upper center", ncol=len(labels), bbox_to_anchor=(0.5, 0.95), frameon=False, fontsize=13)
italicize_species_legend(leg)

plt.tight_layout(rect=(0, 0, 1, 0.92))
save_figure("main/figure2a-c.png")
