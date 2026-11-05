import matplotlib.pyplot as plt
from common import FOCAL_SPECIES, TRAITS, genome_layout, load_results, plot_manhattan, save_figure, set_plot_theme

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
            ax.set_title(f"{trait} in " + r"$\it{" + spp + "}$")
    plt.tight_layout()
    save_figure(path, fig)


manhattan_grid("neg_log10_p", "$-log_{10}$p-value", "manhattan_pvalues.png", show_threshold=True)
manhattan_grid("beta", "beta (effect size)", "manhattan_beta.png", show_threshold=False)
