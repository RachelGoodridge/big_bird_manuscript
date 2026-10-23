import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common import COLOR_MAP, FOCAL_SPECIES, add_bird_info, expand_years_alive, format_year_axis, load_master, save_figure, set_plot_theme, style_legend

set_plot_theme()
TRAITS = {"length": ("Beak Length", "mm"), "depth": ("Beak Depth", "mm"), "width": ("Beak Width", "mm"),
          "size": ("Beak Size (L + D + W)", "mm"), "pointedness": ("Beak Pointedness (L / D)", "ratio")}

MIN_BIRDS = 15  # only show years with at least this many measured birds of a species

# every measured fortis / scandens on Daphne, for every year it was alive
master = load_master()
birds = master[master["Species"].isin(FOCAL_SPECIES)].dropna(subset=["bill.length", "bill.depth", "bill.width"]).rename(columns={"meaningful.unique": "ind"})
birds = add_bird_info(birds, master)
birds = birds.assign(length=birds["bill.length"], depth=birds["bill.depth"], width=birds["bill.width"],
                     size=birds["bill.length"] + birds["bill.depth"] + birds["bill.width"], pointedness=birds["bill.length"] / birds["bill.depth"])
alive = expand_years_alive(birds)
alive = alive[alive.groupby(["Year_alive", "Species"])["ind"].transform("size") >= MIN_BIRDS]

# yearly mean, standard deviation, coefficient of variation ((sd / mean) * 100, as in Grant & Grant), and number of birds
grouped = alive.groupby(["Year_alive", "Species"])[list(TRAITS)]
mean, sd = grouped.mean(), grouped.std()
cv = sd / mean * 100
n_birds = grouped.size()
point_size = 5 + 40 * n_birds / n_birds.max()


############################
# top row: yearly means +/- sd, bottom row: coefficients of variation (point size = number of birds)
fig, axes = plt.subplots(2, len(TRAITS), figsize=(30, 9.75))
for col, (trait, (title, unit)) in enumerate(TRAITS.items()):
    for spp in FOCAL_SPECIES:
        years = mean.xs(spp, level="Species").index
        m, s, c = mean.xs(spp, level="Species")[trait], sd.xs(spp, level="Species")[trait], cv.xs(spp, level="Species")[trait]
        size = point_size.xs(spp, level="Species")
        axes[0, col].errorbar(years, m, yerr=s, color=COLOR_MAP[spp], ecolor=COLOR_MAP[spp], elinewidth=0.8, capsize=2, alpha=0.6, linewidth=0, zorder=1)
        axes[0, col].plot(years, m, color=COLOR_MAP[spp], linewidth=1.5, zorder=2, label=spp)
        axes[0, col].scatter(years, m, s=size, color=COLOR_MAP[spp], edgecolor="none", zorder=3)
        axes[1, col].plot(years, c, color=COLOR_MAP[spp], linewidth=1.5, zorder=2, label=spp)
        axes[1, col].scatter(years, c, s=size, color=COLOR_MAP[spp], edgecolor="none", zorder=3)

    axes[0, col].set_title(title)
    axes[0, col].set_ylabel(f"Mean ({unit})")
    axes[1, col].set_ylabel("Coefficient of Variation (%)")
    axes[1, col].set_ylim(0, cv.max().max() * 1.05)  # same scale for every trait

for ax in axes.flat:
    format_year_axis(ax)
    ax.set_xlim(1982, 2024)
    sns.despine(top=True, right=True, left=False, bottom=False, ax=ax)
style_legend(axes[0, 0], loc="upper right")

plt.tight_layout()
save_figure("supplement/beak_traits_over_time.png", fig)
