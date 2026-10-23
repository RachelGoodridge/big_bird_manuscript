from common import load_allometry_birds, plot_allometry, set_plot_theme

set_plot_theme()

# allometric scaling of beak size with body size, without magnirostris, fuliginosa, hybrids, or unknowns (fit to all birds shown)
birds = load_allometry_birds()
birds = birds[~birds["Species"].isin(["magnirostris", "fuliginosa", "hybrid", "unknown", "Unknown"])]
plot_allometry(birds, birds, r"excluding $\it{magnirostris}$, $\it{fuliginosa}$, and hybrids", "supplement/allometry.png")
