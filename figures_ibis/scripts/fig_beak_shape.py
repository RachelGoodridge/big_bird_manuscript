import matplotlib.pyplot as plt
from common import load_ibd_beak_birds, load_master, plot_beak_vs_ibd, save_figure, set_plot_theme, style_legend

set_plot_theme()

# birds carrying IBD with the founder, with beak measurements
birds = load_ibd_beak_birds(load_master())

# plot fraction of SNPs against beak size and beak shape
fig, axes = plt.subplots(1, 2, figsize=(10.6, 5))
plot_beak_vs_ibd(axes[0], birds, "PC(Beak Size)", show_legend=False)
plot_beak_vs_ibd(axes[1], birds, "PC(Beak Pointedness)", show_legend=True)
style_legend(axes[1], loc="upper right")

plt.tight_layout(rect=[0, 0, 1, 0.88])
save_figure("main/beak_shape_vs_frac_ibd.png")
