import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from common import TRAITS, cell_summary, group_title, load_birds, plot_mean_heatmap, save_figure, set_plot_theme

set_plot_theme()
group, locus1, locus2, trait = "fortis", "G29", "G27", "PC(Beak Size)"

# beak size of fortis in each G29 x G27 genotype combination
birds = load_birds()
birds = birds[birds["Species"] == group].dropna(subset=[locus1, locus2, trait])
means, stds, ns = cell_summary(birds, locus1, locus2, trait, by_ibd=False)

# two-way ANOVA (type II) for an interaction between the two loci
model = smf.ols(f"Q('{trait}') ~ C({locus1}) * C({locus2})", data=birds).fit()
interaction = anova_lm(model, typ=2).loc[f"C({locus1}):C({locus2})"]
print(f"{len(birds)} {group}: interaction F({int(interaction['df'])}, {int(model.df_resid)}) = {interaction['F']:.2f}, p = {interaction['PR(>F)']:.3g}")
title = f"{group_title(group)}  |  {locus1} × {locus2}" + f"\ninteraction p = {interaction['PR(>F)']:.3g}"

plot_mean_heatmap(means, ns, locus1, locus2, trait, title, f"{group}_{locus1}_x_{locus2}_{TRAITS[trait]}_mean.png")

# deviation of each cell mean from a purely additive model (row mean + column mean - grand mean); positive = larger than additive
residuals = means - (means.mean(axis=1).to_numpy()[:, None] + means.mean(axis=0).to_numpy()[None, :] - np.nanmean(means.to_numpy()))

# interaction plot: mean (± SE) beak size across G29 genotypes, one line per G27 genotype, labeled with the residual from additivity
ses = stds / np.sqrt(ns)
fig, ax = plt.subplots(figsize=(6.5, 5))
palette = sns.color_palette("colorblind", n_colors=len(means.columns))
x = np.arange(len(means.index))
for color, g2 in zip(palette, means.columns):
    ax.errorbar(x, means[g2], yerr=ses[g2], label=g2, marker="o", linewidth=2.2, markersize=7, capsize=4, capthick=1.2, color=color)

# residual labels placed clear of the lines: left of the first points and right of the last ones (no lines there; spread apart vertically
# where they would overlap), and below right of the points in between (their lines leave upward to the right)
min_gap = 0.07 * (means.max().max() - means.min().min())
for xi, g1 in enumerate(means.index):
    label_y = means.loc[g1].sort_values()
    for i in range(1, len(label_y)):
        label_y.iloc[i] = max(label_y.iloc[i], label_y.iloc[i - 1] + min_gap)
    if xi == 0:
        place = dict(xytext=(-10, 0), ha="right", va="center")
    elif xi == len(means.index) - 1:
        place = dict(xytext=(10, 0), ha="left", va="center")
    else:
        place = dict(xytext=(9, -5), ha="left", va="top")
    for color, g2 in zip(palette, means.columns):
        y = means.loc[g1, g2] if 0 < xi < len(means.index) - 1 else label_y[g2]
        ax.annotate(f"{residuals.loc[g1, g2]:+.2f}", (xi, y), textcoords="offset points", fontsize="x-small", color=color, **place)
ax.set_xlim(-0.45, len(x) - 0.6)
ax.set_xticks(x, means.index)
ax.xaxis.grid(False)
ax.set_xlabel(locus1)
ax.set_ylabel(trait)
ax.set_title(title)
ax.legend(title=locus2, loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)
sns.despine(top=True, right=True, ax=ax)
plt.tight_layout()
save_figure(f"{group}_{locus1}_x_{locus2}_{TRAITS[trait]}_interaction.png", fig)
