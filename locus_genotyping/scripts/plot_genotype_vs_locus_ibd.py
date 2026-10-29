import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from common import (COLOR_MAP, FOCAL_SPECIES, FOUNDER, GENOTYPES, LOCUS_IBD_COVERAGE, founder_ibd_per_bird, load_genotypes, load_loci, load_master,
                    locus_ibd_fraction, save_figure, set_plot_theme, style_legend)

set_plot_theme()
np.random.seed(0)  # fix the strip plot jitter so the figure is reproducible

# fortis and scandens carrying any IBD with the founder, and the fraction of each locus their segments cover
master = load_master()
birds = founder_ibd_per_bird(master)
birds = birds[birds["Species"].isin(FOCAL_SPECIES) & (birds["fraction_snps_ibd"] > 0)].set_index("ind")
locus_ibd = locus_ibd_fraction(master)

# genotype at each locus vs. the fraction of the locus in IBD; the dashed line is the cutoff for being in IBD at the locus
fig, axes = plt.subplots(2, 3, figsize=(13, 8.5), sharey=True)
for ax, locus in zip(axes.flat, load_loci()["name"]):
    df = birds.assign(genotype=load_genotypes(locus), locus_ibd=locus_ibd[locus]).dropna(subset=["genotype"])
    sns.boxplot(x="genotype", y="locus_ibd", hue="Species", hue_order=FOCAL_SPECIES, palette=COLOR_MAP, data=df, order=GENOTYPES[locus],
                boxprops={"alpha": 0.5}, showfliers=False, legend=False, ax=ax, zorder=1)
    sns.stripplot(x="genotype", y="locus_ibd", hue="Species", hue_order=FOCAL_SPECIES, palette=COLOR_MAP, data=df, order=GENOTYPES[locus],
                  dodge=True, size=5, alpha=0.5, legend=ax is axes.flat[0], ax=ax, zorder=2)
    ax.axhline(LOCUS_IBD_COVERAGE, color="gray", linestyle="--", linewidth=1, zorder=0)
    ax.set_xlabel(None)
    ax.set_ylabel("Fraction of Locus in IBD")
    ax.set_title(f"Locus {locus} (founder {load_genotypes(locus)[FOUNDER]})")
    sns.despine(top=True, right=True, ax=ax)

# one legend for all panels above them
handles, labels = axes.flat[0].get_legend_handles_labels()
axes.flat[0].get_legend().remove()
style_legend(fig, handles=handles, labels=labels, loc="lower center", bbox_to_anchor=(0.5, 0.98), ncol=2, frameon=False, markerscale=1.5)
plt.tight_layout()
save_figure("gts_vs_fracIBD_spp.png", fig)
