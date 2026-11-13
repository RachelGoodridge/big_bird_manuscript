import matplotlib.pyplot as plt
from common import control_yearly_summary, load_master, plot_controls, save_figure, set_plot_theme

set_plot_theme()

# IBD carrier frequency (a) and mean fraction of SNPs in IBD (b) among fortis / scandens over time,
# for the founder vs. each negative control treated as a founder. The labels say plain "IBD", without the subscript F of the other figures,
# since most of these segments are shared with a control bird and not with the founder
summary = control_yearly_summary(load_master())

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12.5, 6))
handles = plot_controls(ax_a, summary, "carrier_freq", "IBD Carrier Frequency")
plot_controls(ax_b, summary, "mean_ibd", "Mean Proportion of SNPs in IBD")

# one legend for both panels above them (columns: founder, mean of controls, each control)
fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.94), frameon=False)
plt.tight_layout(rect=(0, 0, 1, 0.86))
save_figure("supplement/year_vs_IBD_control.png")
