import matplotlib.pyplot as plt
from common import load_allometry_birds, plot_trait_scatter, save_figure, set_plot_theme, style_legend

set_plot_theme()

# body size, beak size and beak pointedness against each other, for the birds of supp_allometry.py (no magnirostris, fuliginosa, hybrids, or
# unknowns), without regressions
PAIRS = [("PC(Body Size)", "PC(Beak Size)"), ("PC(Beak Size)", "PC(Beak Pointedness)"), ("PC(Body Size)", "PC(Beak Pointedness)")]  # x, y of each panel

birds = load_allometry_birds()
birds = birds[~birds["Species"].isin(["magnirostris", "fuliginosa", "hybrid", "unknown", "Unknown"])]

fig, axes = plt.subplots(1, len(PAIRS), figsize=(21, 6))
for ax, (x, y) in zip(axes, PAIRS):
    handles = plot_trait_scatter(ax, birds, x, y)
style_legend(axes[-1], handles=handles, bbox_to_anchor=(1, 1), loc="upper left")
plt.tight_layout()
save_figure("supplement/morphospace.png", fig)
