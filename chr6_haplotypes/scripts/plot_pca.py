import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common import GROUPS, WORK, load_birds, save_figure, set_plot_theme

# PC1 vs PC2 of the LD-pruned SNPs in the chromosome 6 region (genotypes, 02_pca.sh), one color per group; the founder (the same red as the
# conirostris) is drawn last, with a black outline
set_plot_theme()
pca = pd.read_csv(f"{WORK}/pca.eigenvec", sep="\t").rename(columns={"IID": "ind"})[["ind", "PC1", "PC2"]].merge(load_birds(), on="ind")
eigenval = pd.read_csv(f"{WORK}/pca.eigenval", header=None)[0]
pve = 100 * eigenval / eigenval.sum()

fig, ax = plt.subplots(figsize=(8, 6))
for group in ["scandens_no_ibd", "magnirostris_dap", "scandens_ibd", "conirostris_esp", "founder"]:  # the last ones drawn on top
    g = GROUPS.loc[group]
    sub = pca[pca["group"] == group]
    founder = group == "founder"
    ax.scatter(sub["PC1"], sub["PC2"], s=70 if founder else 35, color=g["color"], edgecolor="black" if founder else "none",
               linewidth=1 if founder else 0, label=g["label"])
ax.set_xlabel(f"PC1 ({pve[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({pve[1]:.1f}%)")
ax.set_title("Chromosome 6 (10-30.5 Mb)")
sns.despine(ax=ax)
handles, labels = ax.get_legend_handles_labels()
order = [labels.index(label) for label in GROUPS["label"] if label in labels]  # legend in the groups' order
ax.legend([handles[k] for k in order], [labels[k] for k in order], loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False, markerscale=1)
save_figure("pca.png", fig)

# each group's mean position (where the scandens with IBD sit relative to the others)
print(pca.groupby("group", sort=False)[["PC1", "PC2"]].agg(["count", "mean"]).iloc[:, [0, 1, 3]].round(3).to_string())
