from common import load_allometry_birds, plot_allometry, set_plot_theme

set_plot_theme()

# allometric scaling of beak size with body size for all birds (except unknowns), fit to only the ground finches
birds = load_allometry_birds()
birds = birds[~birds["Species"].isin(["unknown", "Unknown"])]
ground = birds[birds["Species"].isin(["fuliginosa", "fortis", "magnirostris"])]
plot_allometry(birds, ground, r"fit to $\it{fuliginosa}$, $\it{fortis}$, and $\it{magnirostris}$", "supplement/allometry_ground.png")
