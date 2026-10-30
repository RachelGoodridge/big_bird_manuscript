import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from common import GENOTYPING, load_chromosomes, load_drought_fortis, load_master, load_trimmed_segments, save_figure, set_plot_theme, style_legend

set_plot_theme()
birds = load_drought_fortis(load_master())

# trimmed segments with real chromosome names and positions in bp
seg = load_trimmed_segments()
seg["real_chr"] = seg["chrom"].map({code: f"chr{name}" for code, name in load_chromosomes().items()})
seg["phys_start_pos"] *= 1e6
seg["phys_end_pos"] *= 1e6

# the 6 large-effect beak loci
loci_pos = pd.read_csv(f"{GENOTYPING}/data/loci_positions.csv")
loci_pos["name"] = loci_pos["locus"].str.split("_").str[-1]
loci_pos = loci_pos.sort_values("name")

# a bird is in IBD at a locus if its segments cover at least 80% of the locus
ibd_cutoff = 0.8
for _, row in loci_pos.iterrows():
    segs = seg[seg["real_chr"] == row["chrom"]]
    overlap = (segs["phys_end_pos"].clip(upper=row["pos_end"]) - segs["phys_start_pos"].clip(lower=row["pos_start"])).clip(lower=0)
    coverage = overlap.groupby(segs["ind"]).sum() / (row["pos_end"] - row["pos_start"])
    birds[f"ibd_{row['name']}"] = np.where(birds["meaningful.unique"].map(coverage).fillna(0) >= ibd_cutoff, "IBD", "no IBD")

# bar plot of the proportion that died with and without IBD at each locus
# one-sided Fisher's exact test: are individuals in IBD at the locus more likely to have died?
ibd_colors = {"no IBD": "#A9A9A9", "IBD": "#BF5858"}
bar_width = 0.38
fig, ax = plt.subplots(figsize=(9, 4.5))
for i, name in enumerate(loci_pos["name"]):
    counts = pd.crosstab(birds[f"ibd_{name}"], birds["outcome"]).reindex(index=["IBD", "no IBD"], columns=["died", "survived"], fill_value=0)
    _, p = stats.fisher_exact(counts, alternative="greater")
    for j, status in enumerate(["no IBD", "IBD"]):
        n_died, n_total = counts.loc[status, "died"], counts.loc[status].sum()
        x = i + (j - 0.5) * bar_width
        ax.bar(x, n_died / n_total, width=bar_width, color=ibd_colors[status], label=status if i == 0 else None)
        ax.text(x, n_died / n_total + 0.02, f"{n_died}/{n_total}", ha="center", va="bottom", fontsize="x-small")
    ax.text(i, 1.12, f"p = {p:.2g}", ha="center", va="bottom", fontsize="small")
ax.set_xticks(range(len(loci_pos)))
ax.set_xticklabels([f"Locus {name}" for name in loci_pos["name"]])
ax.set_ylim(0, 1.2)
ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.grid(axis="x", visible=False)
ax.set_ylabel("Proportion Died in 2004 Drought")
ax.set_title("One-sided Fisher's exact test (IBD individuals more likely to die)", fontsize="small")
style_legend(ax, loc="upper left", bbox_to_anchor=(1.01, 1))
sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)
plt.tight_layout()
save_figure("barplot_survival_locus_IBD.png")
