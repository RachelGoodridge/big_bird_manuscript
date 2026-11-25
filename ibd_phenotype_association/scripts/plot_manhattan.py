import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from common import FOCAL_SPECIES, TRAITS, genome_layout, load_results, plot_manhattan, save_figure, set_plot_theme, species_label

TICK_NUDGE = {"26": -5, "28": 5}  # horizontal shifts (points) of the x labels of chromosomes 26 and 28, which have too few windows to fit side by side

set_plot_theme()
genome = genome_layout()
df = load_results(genome)


def manhattan_grid(y_col, y_label, path, show_threshold):
    """One panel per trait (rows) and species (columns)."""
    fig, axes = plt.subplots(len(TRAITS), len(FOCAL_SPECIES), figsize=(30, 16))
    for col, spp in enumerate(FOCAL_SPECIES):
        for row, trait in enumerate(TRAITS):
            ax = axes[row, col]
            plot_manhattan(ax, df[(df["species"] == spp) & (df["trait"] == trait)], genome, spp, y_col, y_label, show_threshold)
            ax.set_title(f"{trait} in " + species_label(spp))
    for ax in axes.flat:
        for label in ax.get_xticklabels():
            if label.get_text() in TICK_NUDGE:
                label.set_transform(label.get_transform() + mtransforms.ScaledTranslation(TICK_NUDGE[label.get_text()] / 72, 0, fig.dpi_scale_trans))
    plt.tight_layout()
    save_figure(path, fig)


manhattan_grid("neg_log10_p", "$-log_{10}$p-value", "manhattan_pvalues.png", show_threshold=True)
manhattan_grid("beta", "beta (effect size)", "manhattan_beta.png", show_threshold=False)
