import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from common import COLOR_MAP, IBD_LABEL, load_drought_fortis, load_master, save_figure, set_plot_theme, style_legend

set_plot_theme()
np.random.seed(0)  # fix the strip plot jitter so the figure is reproducible
birds = load_drought_fortis(load_master())

# one-sided Fisher's exact test: are individuals with any IBD more likely to have died?
ibd_counts = pd.crosstab(birds["has_ibd"], birds["outcome"]).reindex(index=[True, False], columns=["died", "survived"], fill_value=0)
_, p_ibd = stats.fisher_exact(ibd_counts, alternative="greater")

# fraction of SNPs in IBD for birds that survived or died in the 2004 drought
fig, ax = plt.subplots(figsize=(6, 5))
sns.boxplot(x="outcome", y="fraction_snps_ibd", data=birds, boxprops={"facecolor": "none"}, ax=ax, showfliers=False)
sns.stripplot(x="outcome", y="fraction_snps_ibd", data=birds, hue="Species", palette=COLOR_MAP, ax=ax)
ax.set_title(f"One-sided Fisher's exact test (IBD individuals more likely to die), p = {p_ibd:.2g}", fontsize=9)
ax.set_xlabel("Survival Outcome in 2004 Drought")
ax.set_ylabel(IBD_LABEL)
style_legend(ax, loc="upper right", bbox_to_anchor=(1.15, 1), title=None)
sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)
plt.tight_layout()
save_figure("supplement/boxplot_survival_IBD.png")
