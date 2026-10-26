import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common import COLOR_MAP, GROUPS, IBD_LABEL, REPO, founder_ibd_per_bird, load_master, order_groups, save_figure, set_plot_theme, style_legend

set_plot_theme()

# fraction of SNPs in IBD with the founder for every bird
birds = order_groups(founder_ibd_per_bird(load_master()), GROUPS)

# conirostris ancestry from ADMIXTURE (K = 9, the conirostris cluster is column 1)
anc = pd.read_csv(f"{REPO}/ancestry_estimates/whole_data_new/admixture/finches.9.Q", header=None, sep="\\s+")
anc["ind"] = pd.read_csv(f"{REPO}/ancestry_estimates/whole_data_new/admixture/finches.fam", header=None, sep="\\s+")[1]
birds["anc_con"] = birds["ind"].map(anc.set_index("ind")[1])

plt.figure(figsize=(7.04, 5.28))
sns.scatterplot(x="anc_con", y="fraction_snps_ibd", hue="Species", palette=COLOR_MAP, data=birds, s=50, alpha=0.7)
plt.xlabel("Conirostris Ancestry Proportion")
plt.ylabel(IBD_LABEL)
sns.despine(top=True, right=True, left=True, bottom=True)
style_legend(plt.gca(), title=None)
plt.tight_layout()
save_figure("supplement/anc_con_vs_frac_ibd.png")
