import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2
from common import FOCAL_SPECIES, RESULTS, TRAITS, load_birds, save_figure, set_plot_theme

set_plot_theme()
df = pd.read_csv(RESULTS, dtype={"chrom": str})
df["neg_log10_p"] = -np.log10(df["p_value"])
panels = [(spp, trait) for spp in FOCAL_SPECIES for trait in TRAITS]


def panel_grid():
    fig, axes = plt.subplots(len(FOCAL_SPECIES), len(TRAITS), figsize=(20, 10))
    return fig, axes.flat


# QQ plots of the observed vs expected (uniform) p-values, with the genomic inflation factor (lambda)
fig, axes = panel_grid()
for ax, (spp, trait) in zip(axes, panels):
    p = df.loc[(df["species"] == spp) & (df["trait"] == trait), "p_value"].to_numpy()
    expected = -np.log10(np.arange(1, len(p) + 1) / len(p))
    observed = -np.log10(np.sort(p))
    inflation = np.median(chi2.isf(p, df=1)) / chi2.ppf(0.5, df=1)
    ax.scatter(expected, observed, s=2, alpha=0.5, color="#1a5276")
    ax.plot([0, expected.max()], [0, expected.max()], color="red", linewidth=1, linestyle="--")
    ax.set_xlabel("Expected $-log_{10}$p")
    ax.set_ylabel("Observed $-log_{10}$p")
    ax.set_title(r"$\it{" + spp + "}$" + f": {trait}, λ = {inflation:.3f}")
    sns.despine(top=True, right=True, ax=ax)
plt.tight_layout()
save_figure("qq_plots.png", fig)

# number of birds in the test group vs the p-value (does power, rather than biology, drive the signal?)
fig, axes = panel_grid()
for ax, (spp, trait) in zip(axes, panels):
    sub = df[(df["species"] == spp) & (df["trait"] == trait)]
    slope, intercept = np.polyfit(sub["n_test"], sub["neg_log10_p"], 1)
    x = np.sort(sub["n_test"].unique())
    ax.scatter(sub["n_test"], sub["neg_log10_p"], s=2, alpha=0.5, color="#1a5276")
    ax.plot(x, slope * x + intercept, color="red", linewidth=1)
    ax.set_xlabel("Birds in the test (IBD) group")
    ax.set_ylabel("$-log_{10}$p-value")
    ax.set_title(r"$\it{" + spp + "}$" + f": {trait}")
    sns.despine(top=True, right=True, ax=ax)
plt.tight_layout()
save_figure("n_test_vs_pvalue.png", fig)

# sign check of the PCs fit within each species: each should increase with the measurement it summarizes
checks = [("PC(Body Size)", lambda b: b["weight"], "Weight"),
          ("PC(Beak Size)", lambda b: b[["bill.length", "bill.depth", "bill.width"]].sum(axis=1), "Beak length + depth + width"),
          ("PC(Beak Pointedness)", lambda b: b["bill.length"] / b["bill.depth"], "Beak length / depth")]
fig, axes = plt.subplots(len(FOCAL_SPECIES), len(checks), figsize=(12, 4 * len(FOCAL_SPECIES)))
for row, spp in enumerate(FOCAL_SPECIES):
    birds = load_birds(spp)
    for ax, (pc, measure, x_label) in zip(axes[row], checks):
        x = measure(birds)
        ax.scatter(x, birds[pc], s=3, alpha=0.3, color="#1a5276")
        ax.set_xlabel(x_label)
        ax.set_ylabel(pc)
        ax.set_title(r"$\it{" + spp + "}$" + f": r = {np.corrcoef(x, birds[pc])[0, 1]:.3f}")
        sns.despine(top=True, right=True, ax=ax)
plt.tight_layout()
save_figure("pc_sign_check.png", fig)
