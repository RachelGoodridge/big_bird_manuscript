import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from common import COLOR_MAP, load_drought_fortis, load_master, save_figure, set_plot_theme, style_legend

set_plot_theme()
np.random.seed(0)  # fix the strip plot jitter so the figure is reproducible
birds = load_drought_fortis(load_master())

# Mann-Whitney U test on beak size between outcomes
_, p_beak = stats.mannwhitneyu(*[g["beak_sum"].dropna() for _, g in birds.groupby("outcome")])

# beak size for birds that survived or died in the 2004 drought
fig, ax = plt.subplots(figsize=(6, 5))
sns.boxplot(x="outcome", y="beak_sum", data=birds, boxprops={"facecolor": "none"}, ax=ax, showfliers=False)
sns.stripplot(x="outcome", y="beak_sum", data=birds, hue="Species", palette=COLOR_MAP, ax=ax)
ax.set_title(f"Mann-Whitney U test, p = {p_beak:.2g}", fontsize=9)
ax.set_xlabel("Survival Outcome in 2004 Drought")
ax.set_ylabel("Beak size (length + depth + width)")
style_legend(ax, loc="upper right", bbox_to_anchor=(1.15, 1), title=None)
sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)
plt.tight_layout()
save_figure("supplement/boxplot_survival_beak_size.png")
