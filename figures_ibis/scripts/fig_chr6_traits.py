import importlib.util
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from common import ASSOCIATION, COLOR_MAP, load_chromosomes, save_figure, set_plot_theme

# the IBD-phenotype association helpers (loaded by path since it is also named common.py)
_spec = importlib.util.spec_from_file_location("association_common", f"{ASSOCIATION}/scripts/common.py")
sys.modules["association_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["association_common"])
from association_common import TRAITS, chromosome_lmm, chromosome_test_birds

set_plot_theme()
spp = "scandens"
chrom_name = "6"

# ── birds + IBD ───────────────────────────────────────────────────────────────
# the scandens of the IBD-phenotype association test (all six measurements and known sex, traits fit within them), with the fraction of
# chromosome 6 SNPs in IBD (chr_ibd) and of the SNPs outside it (loco_ibd), and the chromosome 6 LOCO GRM
chrom_code = {name: code for code, name in load_chromosomes().items()}[chrom_name]
birds, eig = chromosome_test_birds(spp, chrom_code)

# create IBD categories for boxplot
thresh = 0.5
IBD_type = lambda x: f"IBD > {thresh}" if x > thresh else ("Low IBD" if x > 0 else "No IBD")
birds["IBD_type"] = birds["chr_ibd"].apply(IBD_type)
birds_box = birds[birds["IBD_type"] != "Low IBD"].copy()

# ── regression stats (GRM-LMM on the continuous chr_ibd before binning) ──────
# the same LOCO GRM, covariates (sex, IBD outside the chromosome, and PC(Body Size) for beak size / pointedness), and likelihood ratio test
# as the windows of the IBD-phenotype association, but with continuous IBD over the whole chromosome
stats = {trait: dict(zip(["beta", "p_reg"], chromosome_lmm(birds, eig, trait))) for trait in TRAITS}

# ── Mann-Whitney (boxplot groups) ─────────────────────────────────────────────
for trait in TRAITS:
    a = birds_box.loc[birds_box["IBD_type"] == "No IBD", trait]
    b = birds_box.loc[birds_box["IBD_type"] == f"IBD > {thresh}", trait]
    _, p_mwu = mannwhitneyu(a, b, alternative="two-sided")
    stats[trait]["p_mwu"] = p_mwu
print(birds_box["IBD_type"].value_counts().to_dict())
print(pd.DataFrame(stats).T)

# ── plot ──────────────────────────────────────────────────────────────────────
def annotate_panel(ax, beta, p_reg, p_mwu, y_max, y_range):
    # bracket: Mann-Whitney stars
    if p_mwu < 0.001:   stars = "***"
    elif p_mwu < 0.01:  stars = "**"
    elif p_mwu < 0.05:  stars = "*"
    else:               stars = "ns"

    h = y_range * 0.04
    pad = y_range * 0.02
    y_top = y_max + y_range * 0.05

    ax.plot([0, 0, 1, 1], [y_top, y_top + h, y_top + h, y_top], lw=1.2, color="black", clip_on=False)
    ax.text(0.5, y_top + h + pad, stars, ha="center", va="bottom", fontsize=18)
    ax.set_ylim(top=y_top + h + pad + y_range * 0.08)

    # regression stats (continuous chr_ibd, GRM-LMM) placed below the panel - this is a
    # separate calculation from the boxed/binned Mann-Whitney comparison shown above
    p_reg_str = f"p = {p_reg:.3f}" if p_reg > 0.001 else "p < 0.001"
    ax.text(0.5, -0.14, f"GRM-LMM: β = {beta:.3f}, {p_reg_str}", transform=ax.transAxes, ha="center", va="top", fontsize="small", color="#333333")

fig, axes = plt.subplots(2, 2, figsize=(10.6, 11))
for ax, trait in zip(axes.flat, TRAITS):
    sns.boxplot(data=birds_box, x="IBD_type", y=trait, order=["No IBD", f"IBD > {thresh}"], color=COLOR_MAP[spp], zorder=2, ax=ax)
    if ax.collections:
        ax.collections[-1].set_zorder(3)

    ax.set_xlabel(None)
    ax.set_ylabel(trait)
    ax.set_title(f"IBD on chr. {chrom_name} in " + r"$\it{" + spp + "}$")

    y_max = birds_box[trait].max()
    y_range = birds_box[trait].max() - birds_box[trait].min()
    annotate_panel(ax, stats[trait]["beta"], stats[trait]["p_reg"], stats[trait]["p_mwu"], y_max, y_range)
    sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)

plt.tight_layout()
plt.subplots_adjust(hspace=0.5)  # room for the regression annotation below each panel
save_figure(f"main/chr{chrom_name}_IBD_vs_traits.png")
