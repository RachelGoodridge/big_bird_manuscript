import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from common import (COLOR_MAP, add_beak_body_ratio, add_bill_pcs, add_body_pc, eigen_grm, lmm_lrt, load_bim, load_chromosomes, load_founder_coef,
                    load_loco_grm, load_master, load_trimmed_segments, save_figure, set_plot_theme)

set_plot_theme()
spp = "scandens"
chrom_name = "6"

# ── SNP counts ────────────────────────────────────────────────────────────────
# number of SNPs on the chromosome of interest vs. the rest of the genome
chrom_code = {name: code for code, name in load_chromosomes().items()}[chrom_name]
bim = load_bim()
chr_snps = (bim["chrom"] == chrom_code).sum()
background_snps = len(bim) - chr_snps

# ── metadata ──────────────────────────────────────────────────────────────────
# birds in the IBIS run (sex coded male = 0, female = 1)
master = load_master()
ibis_birds = set(load_founder_coef()["ind"])
master2 = master[master["meaningful.unique"].isin(ibis_birds) & (master["Species"] == spp)].copy()

# ── IBD: chromosome of interest + background (chromosome excluded) ────────────
seg = load_trimmed_segments()
chr_IBD = seg[seg["chrom"] == chrom_code].groupby("ind")["SNP_count"].sum()
background_IBD = seg[seg["chrom"] != chrom_code].groupby("ind")["SNP_count"].sum()

# fraction of the SNPs on the chromosome / in the background that are in IBD (same measure as the other figures)
master2["chr_IBD"] = master2["meaningful.unique"].map(chr_IBD).fillna(0) / chr_snps
master2["background_IBD"] = master2["meaningful.unique"].map(background_IBD).fillna(0) / background_snps

# ── traits + IBD category ─────────────────────────────────────────────────────
# remove NA values for phenotypes
master2 = master2.dropna(subset=["weight", "wing", "tarsus", "bill.length", "bill.depth", "bill.width"])

# PCA on body size (signed so larger values mean heavier birds) and on bill dimensions
master2 = add_bill_pcs(add_body_pc(master2))

# beak volume relative to body mass, standardized within these scandens
master2 = add_beak_body_ratio(master2)

# create IBD categories for boxplot
thresh = 0.5
IBD_type = lambda x: f"IBD > {thresh}" if x > thresh else ("Low IBD" if x > 0 else "No IBD")
master2["IBD_type"] = master2["chr_IBD"].apply(IBD_type)
master2_box = master2[master2["IBD_type"] != "Low IBD"].copy()

# ── regression stats (GRM-LMM on full continuous chr_IBD before binning) ─────
# uses the same LOCO GRM and LMM as the IBD-phenotype association (ibd_phenotype_association/), so this
# is a consistent statistical zoom-in on the same signal shown in the Manhattan plot
trait_list = ["PC(Body Size)", "PC(Beak Size)", "PC(Beak Pointedness)", "Beak:Body Ratio"]

# beak size/pointedness additionally control for body size (allometry); background_IBD
# also LOCO is always included as a fixed-effect covariate
BASE_COVARS = ["sex", "background_IBD"]
EXTRA_COVARS = {
    "PC(Beak Size)": ["PC(Body Size)"],
    "PC(Beak Pointedness)": ["PC(Body Size)"],
}

grm_matrix, grm_ids = load_loco_grm(spp, chrom_code)
grm_index = {iid: i for i, iid in enumerate(grm_ids)}

# regression uses a subset restricted to individuals present in the LOCO GRM -
# the box plots below are untouched and keep the full master2 sample
master2_grm = master2[master2["meaningful.unique"].isin(grm_index)].copy()
missing_from_grm = set(master2["meaningful.unique"]) - set(grm_index)
if missing_from_grm:
    print(f"[warning] {len(missing_from_grm)}/{len(master2)} {spp} individuals not in the "
          f"LOCO GRM and excluded from the regression: {sorted(missing_from_grm)}")

stats = {}
grm_idx = [grm_index[iid] for iid in master2_grm["meaningful.unique"].tolist()]
eig = eigen_grm(grm_matrix[np.ix_(grm_idx, grm_idx)])

for trait in trait_list:
    covar_cols = BASE_COVARS + EXTRA_COVARS.get(trait, [])
    X = np.column_stack([np.ones(len(master2_grm)), master2_grm["chr_IBD"].values] + [master2_grm[c].values for c in covar_cols])
    y = master2_grm[trait].values
    beta, p_reg = lmm_lrt(y, X, eig)
    stats[trait] = {"beta": beta, "p_reg": p_reg}

# ── Mann-Whitney (boxplot groups) ─────────────────────────────────────────────
for trait in trait_list:
    a = master2_box.loc[master2_box["IBD_type"] == "No IBD", trait]
    b = master2_box.loc[master2_box["IBD_type"] == f"IBD > {thresh}", trait]
    _, p_mwu = mannwhitneyu(a, b, alternative="two-sided")
    stats[trait]["p_mwu"] = p_mwu
print(master2_box["IBD_type"].value_counts().to_dict())
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

    # regression stats (continuous chr_IBD, GRM-LMM) placed below the panel - this is a
    # separate calculation from the boxed/binned Mann-Whitney comparison shown above
    p_reg_str = f"p = {p_reg:.3f}" if p_reg > 0.001 else "p < 0.001"
    ax.text(0.5, -0.14, f"GRM-LMM: β = {beta:.3f}, {p_reg_str}", transform=ax.transAxes, ha="center", va="top", fontsize="small", color="#333333")

fig, axes = plt.subplots(2, 2, figsize=(10.6, 11))
for ax, trait in zip(axes.flat, trait_list):
    sns.boxplot(data=master2_box, x="IBD_type", y=trait, order=["No IBD", f"IBD > {thresh}"], color=COLOR_MAP[spp], zorder=2, ax=ax)
    if ax.collections:
        ax.collections[-1].set_zorder(3)

    ax.set_xlabel(None)
    ax.set_ylabel(trait)
    ax.set_title(f"IBD on chr. {chrom_name} in " + r"$\it{" + spp + "}$")

    y_max = master2_box[trait].max()
    y_range = master2_box[trait].max() - master2_box[trait].min()
    annotate_panel(ax, stats[trait]["beta"], stats[trait]["p_reg"], stats[trait]["p_mwu"], y_max, y_range)
    sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)

plt.tight_layout()
plt.subplots_adjust(hspace=0.5)  # room for the regression annotation below each panel
save_figure(f"main/chr{chrom_name}_IBD_vs_traits.png")
