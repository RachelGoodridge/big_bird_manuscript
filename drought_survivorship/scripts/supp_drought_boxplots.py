import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from common import COLOR_MAP, IBD_LABEL, load_drought_fortis, load_master, save_figure, set_plot_theme

set_plot_theme()
np.random.seed(0)  # fix the strip plot jitter so the figure is reproducible
birds = load_drought_fortis(load_master())

# (a) Mann-Whitney U test on beak size between outcomes
_, p_beak = stats.mannwhitneyu(*[g["beak_sum"].dropna() for _, g in birds.groupby("outcome")])

# (b) one-sided Fisher's exact test: are individuals with any IBD more likely to have died?
ibd_counts = pd.crosstab(birds["has_ibd"], birds["outcome"]).reindex(index=[True, False], columns=["died", "survived"], fill_value=0)
_, p_ibd = stats.fisher_exact(ibd_counts, alternative="greater")

# beak size (a) and fraction of SNPs in IBD (b) for birds that survived or died in the 2004 drought
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 6))
for ax, y_col, y_label, title in [
        (ax_a, "beak_sum", "Beak size (length + depth + width)", f"Mann-Whitney U test, p = {p_beak:.2g}"),
        (ax_b, "fraction_snps_ibd", IBD_LABEL, f"One-sided Fisher's exact test\n(IBD individuals more likely to die), p = {p_ibd:.2g}")]:
    sns.boxplot(x="outcome", y=y_col, data=birds, boxprops={"facecolor": "none"}, ax=ax, showfliers=False)
    sns.stripplot(x="outcome", y=y_col, data=birds, hue="Species", palette=COLOR_MAP, ax=ax)
    ax.set_title(title, fontsize="small")
    ax.set_xlabel("Survival Outcome in 2004 Drought")
    ax.set_ylabel(y_label)
    ax.get_legend().remove()
    sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)

# no legend: every bird in the drought data is fortis
plt.tight_layout()
save_figure("supplement/boxplot_survival_beak_size_IBD.png")
