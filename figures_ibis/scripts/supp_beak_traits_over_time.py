import matplotlib.pyplot as plt
import seaborn as sns
from common import (COLOR_MAP, FOCAL_SPECIES, format_year_axis, load_master, load_measured_fortis_scandens, plot_trait_means, save_figure,
                    set_plot_theme, style_legend, yearly_trait_summary)

set_plot_theme()
TRAITS = {"length": ("Beak Length", "mm"), "depth": ("Beak Depth", "mm"), "width": ("Beak Width", "mm"),
          "size": ("Beak Size (L + D + W)", "mm"), "pointedness": ("Beak Pointedness (L / D)", "ratio")}

# every measured fortis / scandens on Daphne, for every year it was alive
birds = load_measured_fortis_scandens(load_master())
birds = birds.assign(length=birds["bill.length"], depth=birds["bill.depth"], width=birds["bill.width"],
                     size=birds["bill.length"] + birds["bill.depth"] + birds["bill.width"], pointedness=birds["bill.length"] / birds["bill.depth"])

# yearly mean, standard deviation, coefficient of variation, and point size for the number of birds (years with at least 15 birds of a species)
mean, sd, cv, point_size = yearly_trait_summary(birds, TRAITS)


############################
# top row: yearly means +/- sd, bottom row: coefficients of variation (point size = number of birds)
fig, axes = plt.subplots(2, len(TRAITS), figsize=(30, 9.75))
for col, (trait, (title, unit)) in enumerate(TRAITS.items()):
    plot_trait_means(axes[0, col], mean, sd, point_size, trait)
    for spp in FOCAL_SPECIES:
        c = cv.xs(spp, level="Species")[trait]
        axes[1, col].plot(c.index, c, color=COLOR_MAP[spp], linewidth=1.5, zorder=2, label=spp)
        axes[1, col].scatter(c.index, c, s=point_size.xs(spp, level="Species"), color=COLOR_MAP[spp], edgecolor="none", zorder=3)

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
