import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
import seaborn as sns
from common import COLOR_MAP, FOCAL_SPECIES, add_bird_info, expand_years_alive, load_master, load_trimmed_segments, save_figure, set_plot_theme, style_legend

set_plot_theme()

# trimmed segments of fortis and scandens for every year each bird was alive
seg = expand_years_alive(add_bird_info(load_trimmed_segments(), load_master()))
seg = seg[seg["Species"].isin(FOCAL_SPECIES)]
years = sorted(seg["Year_alive"].unique())

# skinny box plots of segment length per year, side by side for each species, with outliers as points
fig, ax = plt.subplots(figsize=(9.5, 6))
for j, spp in enumerate(FOCAL_SPECIES):
    color = COLOR_MAP[spp]
    for i, year in enumerate(years):
        lengths = seg[(seg["Year_alive"] == year) & (seg["Species"] == spp)]["length_Mb"].to_numpy()
        if len(lengths) == 0:
            continue
        x = i + (j - 0.5) * 0.3
        ax.boxplot(lengths, positions=[x], widths=0.25, patch_artist=True, showfliers=False, showcaps=False,
                   boxprops={"facecolor": to_rgba(color, 0.4), "edgecolor": color, "linewidth": 0.85},
                   whiskerprops={"color": color, "linewidth": 0.85}, medianprops={"color": color, "linewidth": 0.85})
        q1, q3 = np.percentile(lengths, [25, 75])
        outliers = lengths[(lengths < q1 - 1.5 * (q3 - q1)) | (lengths > q3 + 1.5 * (q3 - q1))]
        ax.scatter(np.full(len(outliers), x), outliers, s=4, color=color, alpha=0.6, linewidth=0)

ax.set_xticks(range(len(years)))
ax.set_xticklabels(years, rotation=90)
ax.set_xlim(-0.6, len(years) - 0.4)
ax.set_ylabel("IBD Segment Length (Mb)")
ax.grid(axis="x", visible=False)
style_legend(ax, handles=[Patch(facecolor=to_rgba(COLOR_MAP[s], 0.4), edgecolor=COLOR_MAP[s], label=s) for s in FOCAL_SPECIES],
             loc="upper right")
sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)

plt.tight_layout()
save_figure("supplement/IBD_len_boxplot_over_time.png")
