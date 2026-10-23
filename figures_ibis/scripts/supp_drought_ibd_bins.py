import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from common import load_drought_fortis, load_master, save_figure, set_plot_theme

set_plot_theme()
birds = load_drought_fortis(load_master())

# proportion that died in the 2004 drought with and without genome-wide IBD
bin_labels = ["No IBD", "IBD"]
birds["ibd_bin"] = pd.cut(birds["fraction_snps_ibd"], bins=[-np.inf, 0, np.inf], labels=bin_labels)
bin_counts = pd.crosstab(birds["ibd_bin"], birds["outcome"]).reindex(index=bin_labels, columns=["survived", "died"], fill_value=0)

# one-sided Fisher's exact test: are individuals with IBD more likely to have died?
_, p_bins = stats.fisher_exact(bin_counts.loc[["IBD", "No IBD"], ["died", "survived"]], alternative="greater")

fig, ax = plt.subplots(figsize=(5, 4.5))
for i, (label, color) in enumerate(zip(bin_labels, ["#A9A9A9", "#BF5858"])):
    n_died, n_total = bin_counts.loc[label, "died"], bin_counts.loc[label].sum()
    ax.bar(i, n_died / n_total, width=0.7, color=color)
    ax.text(i, n_died / n_total + 0.02, f"{n_died}/{n_total}", ha="center", va="bottom", fontsize=9)
ax.set_xticks(range(len(bin_labels)))
ax.set_xticklabels(bin_labels)
ax.set_ylim(0, 1.1)
ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.grid(axis="x", visible=False)
ax.set_ylabel("Proportion Died in 2004 Drought")
ax.set_title(f"One-sided Fisher's exact test, p = {p_bins:.2g}", fontsize=9)
sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)
plt.tight_layout()
save_figure("supplement/barplot_survival_IBD_bins.png")
