import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from common import COLOR_MAP, IBD_LABEL, degree_cutoffs, founder_ibd_per_bird, load_all_coef, load_master, save_figure, set_plot_theme, style_legend

set_plot_theme()
lineages = ["A line", "B line", "Big Birds", "founder"]

# the 20 scandens first seen in 2022 sharing the most IBD with the founder
master = load_master()
birds = founder_ibd_per_bird(master)
birds = birds[birds["fraction_snps_ibd"] > 0].sort_values(by="ind")
top_scandens = birds[(birds["Year"] == 2022) & (birds["Species"] == "scandens")].sort_values("fraction_snps_ibd", ascending=False).head(20)

# kinship between all pairs up to 4th degree, with the species of both birds
coef = load_all_coef()
for i in [1, 2]:
    coef[f"Species{i}"] = coef[f"Individual{i}"].map(master.set_index("meaningful.unique")["Species"])
subset = coef[coef["Degree"].isin([0, 1, 2, 3, 4])].copy()

# cutoffs between 2nd/3rd and 3rd/4th degree relatives
cutoffs = degree_cutoffs(coef, [2, 3])

plt.figure(figsize=(6.4, 4.8))
for _, bird in top_scandens.iterrows():
    # kinship of this bird with the founder and his descendants
    rel = subset[(subset["Individual1"] == bird["ind"]) | (subset["Individual2"] == bird["ind"])].copy()
    rel["spp"] = np.where(rel["Individual1"] == bird["ind"], rel["Species2"], rel["Species1"])
    rel = rel[rel["spp"].isin(lineages)].copy()
    rel["fraction_snps_ibd"] = bird["fraction_snps_ibd"]
    sns.scatterplot(x="fraction_snps_ibd", y="Kinship_Coefficient", hue="spp", data=rel[rel["spp"] != "founder"], palette=COLOR_MAP, s=25, alpha=0.5, legend=False, zorder=1)
    sns.scatterplot(x="fraction_snps_ibd", y="Kinship_Coefficient", data=rel[rel["spp"] == "founder"], color=COLOR_MAP["founder"], s=25, alpha=0.75, legend=False, zorder=2)

for c in cutoffs:
    plt.axhline(y=c, color="gray", linestyle="--", linewidth=1)
plt.xlabel(IBD_LABEL)
plt.ylabel("Kinship coefficient with A/B line or founder")
handles = [Line2D([], [], marker="o", color="w", markerfacecolor=COLOR_MAP[label], markersize=5, label=label) for label in lineages]
style_legend(plt.gca(), handles=handles, loc="upper left")
sns.despine(top=True, right=True, left=False, bottom=False)
plt.tight_layout()
save_figure("supplement/frac_ibd_vs_kinship_top20scandens2022.png")
