import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import linalg
from scipy.optimize import minimize_scalar
from scipy.stats import mannwhitneyu, chi2
from sklearn.preprocessing import StandardScaler
from common import COLOR_MAP, add_bill_pcs, add_body_pc, load_bim, load_chromosomes, load_founder_coef, load_master, load_trimmed_segments, save_figure, set_plot_theme

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
# birds in the IBIS run, with consensus sex (Male = 0, Female = 1)
master = load_master()
master["sex.consensus"] = master["sex.consensus"].fillna(master["Sex"].map({1: "Male", 2: "Female"}))
master["Sex"] = master["sex.consensus"].map({"Male": 0, "Female": 1})
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

# beak volume proxy : body mass, log transform to create normal distribution, and then standardize to remove the unit problem
beak_volume = master2["bill.length"] * master2["bill.depth"] * master2["bill.width"]
master2["Beak:Body Ratio"] = np.log(beak_volume / master2["weight"])
master2["Beak:Body Ratio"] = StandardScaler().fit_transform(master2[["Beak:Body Ratio"]]).ravel()

# create IBD categories for boxplot
thresh = 0.5
IBD_type = lambda x: f"IBD > {thresh}" if x > thresh else ("Low IBD" if x > 0 else "No IBD")
master2["IBD_type"] = master2["chr_IBD"].apply(IBD_type)
master2_box = master2[master2["IBD_type"] != "Low IBD"].copy()

# ── regression stats (GRM-LMM on full continuous chr_IBD before binning) ─────
# uses the same LOCO GRM and LMM framework as QTL_GWAS.py, so this
# is a consistent statistical zoom-in on the same signal shown in the Manhattan plot
trait_list = ["PC(Body Size)", "PC(Beak Size)", "PC(Beak Pointedness)", "Beak:Body Ratio"]

# beak size/pointedness additionally control for body size (allometry); background_IBD
# also LOCO is always included as a fixed-effect covariate
BASE_COVARS = ["Sex", "background_IBD"]
EXTRA_COVARS = {
    "PC(Beak Size)": ["PC(Body Size)"],
    "PC(Beak Pointedness)": ["PC(Body Size)"],
}

grm_dir = f"/SSD/reg259/IBIS/allele_freqs/QTL_GWAS/GRM/{spp}"
grm_fam = pd.read_csv(f"{grm_dir}/plink.fam", header=None, sep=r"\s+", names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"], dtype=str)
grm_index = {iid: i for i, iid in enumerate(grm_fam["IID"].tolist())}
grm_matrix = np.loadtxt(f"{grm_dir}/gemma_grm_excl_chr{chrom_code}.cXX.txt")

# regression uses a subset restricted to individuals present in the LOCO GRM -
# the box plots below are untouched and keep the full master2 sample
master2_grm = master2[master2["meaningful.unique"].isin(grm_index)].copy()
missing_from_grm = set(master2["meaningful.unique"]) - set(grm_index)
if missing_from_grm:
    print(f"[warning] {len(missing_from_grm)}/{len(master2)} {spp} individuals not in the "
          f"LOCO GRM and excluded from the regression: {sorted(missing_from_grm)}")

n_missing_sex = master2_grm["Sex"].isna().sum()
if n_missing_sex:
    print(f"[warning] {n_missing_sex} individuals missing Sex and excluded from the regression")
master2_grm = master2_grm.dropna(subset=["Sex"])


def fit_lmm(y, X, K):
    n = len(y)
    null_cols = [0] + list(range(2, X.shape[1]))

    eigenvalues, eigenvectors = np.linalg.eigh(K)
    Qty = eigenvectors.T @ y
    QtX = eigenvectors.T @ X

    def log_likelihood(delta):
        d = eigenvalues + delta
        dinv = 1.0 / d
        XtHinvX = (QtX * dinv[:, None]).T @ QtX
        XtHinvy = (QtX * dinv[:, None]).T @ Qty
        try:
            beta = linalg.solve(XtHinvX, XtHinvy, assume_a="pos")
        except linalg.LinAlgError:
            return np.inf
        resid = Qty - QtX @ beta
        sigma2 = np.sum(resid**2 * dinv) / n
        logdet = np.sum(np.log(d))
        ll = -0.5 * (n * np.log(2 * np.pi * sigma2) + logdet + n)
        return -ll

    result = minimize_scalar(log_likelihood, bounds=(1e-4, 1e4), method="bounded")
    d_opt = eigenvalues + result.x
    dinv_opt = 1.0 / d_opt

    def get_beta_ll(QtX_use):
        XtHinvX = (QtX_use * dinv_opt[:, None]).T @ QtX_use
        XtHinvy = (QtX_use * dinv_opt[:, None]).T @ Qty
        beta = linalg.solve(XtHinvX, XtHinvy, assume_a="pos")
        resid = Qty - QtX_use @ beta
        sigma2 = np.sum(resid**2 * dinv_opt) / n
        logdet = np.sum(np.log(d_opt))
        ll = -0.5 * (n * np.log(2 * np.pi * sigma2) + logdet + n)
        return beta, ll

    beta_full, ll_full = get_beta_ll(QtX)
    _, ll_null = get_beta_ll(QtX[:, null_cols])
    lrt = max(2 * (ll_full - ll_null), 0)
    return beta_full[1], chi2.sf(lrt, df=1)


stats = {}
grm_idx = [grm_index[iid] for iid in master2_grm["meaningful.unique"].tolist()]
K_sub = grm_matrix[np.ix_(grm_idx, grm_idx)]

for trait in trait_list:
    covar_cols = BASE_COVARS + EXTRA_COVARS.get(trait, [])
    X = np.column_stack([np.ones(len(master2_grm)), master2_grm["chr_IBD"].values] + [master2_grm[c].values for c in covar_cols])
    y = master2_grm[trait].values
    beta, p_reg = fit_lmm(y, X, K_sub)
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
    ax.text(0.5, y_top + h + pad, stars, ha="center", va="bottom", fontsize=14)
    ax.set_ylim(top=y_top + h + pad + y_range * 0.08)

    # regression stats (continuous chr_IBD, GRM-LMM) placed below the panel - this is a
    # separate calculation from the boxed/binned Mann-Whitney comparison shown above
    p_reg_str = f"p = {p_reg:.3f}" if p_reg > 0.001 else "p < 0.001"
    ax.text(0.5, -0.14, f"GRM-LMM: β = {beta:.3f}, {p_reg_str}", transform=ax.transAxes, ha="center", va="top", fontsize=9, color="#333333")

fig, axes = plt.subplots(2, 2, figsize=(10.6, 11))
for ax, trait in zip(axes.flat, trait_list):
    sns.boxplot(data=master2_box, x="IBD_type", y=trait, order=["No IBD", f"IBD > {thresh}"], color=COLOR_MAP[spp], zorder=2, ax=ax)
    if ax.collections:
        ax.collections[-1].set_zorder(3)

    ax.set_xlabel(None)
    ax.tick_params(axis="x", labelsize=12)
    ax.set_ylabel(trait, fontsize=12)
    ax.set_title(f"IBD on chr. {chrom_name} in " + r"$\it{" + spp + "}$", fontsize=12)

    y_max = master2_box[trait].max()
    y_range = master2_box[trait].max() - master2_box[trait].min()
    annotate_panel(ax, stats[trait]["beta"], stats[trait]["p_reg"], stats[trait]["p_mwu"], y_max, y_range)
    sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)

plt.tight_layout()
plt.subplots_adjust(hspace=0.5)  # room for the regression annotation below each panel
save_figure(f"main/chr{chrom_name}_IBD_vs_traits.png")
