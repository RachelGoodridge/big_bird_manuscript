import matplotlib.pyplot as plt
from common import HEATMAP_SNPS, plot_peak_heatmap, save_figure, set_plot_theme, species_label
from ibis_common import ibd_label

# heatmap of every bird's genotype (copies of the donor allele) at the SNPs in the FST peak, with a tree of the birds at the side; the heatmap
# itself (and its settings) is in common.py, as it is also a panel of figures_ibis/scripts/fig_4_combined.py
set_plot_theme()
fig, (ax_t, ax_g, ax) = plt.subplots(1, 3, figsize=(15.5, 13), gridspec_kw={"width_ratios": [0.12, 0.025, 1], "wspace": 0.01})
n_snps, n_shown, n_scandens = plot_peak_heatmap(ax_t, ax_g, ax)
ax.set_title(f"Genotypes in the chr. 6 FST peak ({n_snps} {HEATMAP_SNPS} SNPs; {n_shown} of the {n_scandens} {species_label('scandens')} without "
             f"{ibd_label('IBD')} shown)", pad=28)
save_figure("peak_heatmap.png", fig)
