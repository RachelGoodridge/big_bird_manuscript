import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from sklearn.mixture import GaussianMixture
from common import DATA, GENOTYPE_CONFIDENCE, GENOTYPES, GENOTYPING, OUTPUT, load_loci, load_master, save_figure, set_plot_theme

set_plot_theme()
os.makedirs(OUTPUT, exist_ok=True)
GROUPS = ["fortis", "scandens", "Big Birds", "A line", "B line", "founder"]
SPECIES_MARKERS = {"fortis": "o", "scandens": "s", "Big Birds": "D", "A line": "v", "B line": "^", "founder": "X"}
GENOTYPE_COLORS = {"LL": "tab:blue", "SL": "tab:orange", "SS": "tab:green", "BB": "tab:blue", "PB": "tab:orange", "PP": "tab:green"}

master = load_master()
published = pd.read_csv(f"{DATA}/science.adf6218_table_s1.csv").set_index("sample")


def classify_with_gmm(df, classes, features=("PC1", "PC2")):
    """Fit one Gaussian per published genotype on PC1/PC2, then give every bird without a published genotype the most likely one."""
    known = df["genotype"].notna()
    gmms = {g: GaussianMixture(n_components=1, covariance_type="full" if (df["genotype"] == g).sum() >= 3 else "diag")
               .fit(df.loc[df["genotype"] == g, list(features)]) for g in classes}
    log_probs = np.column_stack([gmms[g].score_samples(df.loc[~known, list(features)]) for g in classes])
    probs = np.exp(log_probs - log_probs.max(axis=1, keepdims=True))
    probs /= probs.sum(axis=1, keepdims=True)

    # posterior of each genotype for the unknowns; published genotypes are kept with confidence 1
    for i, g in enumerate(classes):
        df.loc[~known, f"prob_{g}"] = probs[:, i]
    df["predicted_genotype"] = df["genotype"]
    df.loc[~known, "predicted_genotype"] = np.array(classes)[probs.argmax(axis=1)]
    df["confidence"] = 1.0
    df.loc[~known, "confidence"] = probs.max(axis=1)
    return df


for locus in load_loci()["name"]:
    # PCA of the birds' genotypes at the locus SNPs, for Daphne fortis, scandens, and the founder's lineages
    pca = pd.read_csv(f"{GENOTYPING}/work/pca_{locus}.eigenvec", sep=r"\s+").rename(columns={"#IID": "IID"})
    eigenval = pd.read_csv(f"{GENOTYPING}/work/pca_{locus}.eigenval", header=None)[0]
    pca["Species"] = pca["IID"].map(master.set_index("meaningful.unique")["Species"])
    pca = pca[pca["Species"].isin(GROUPS)].reset_index(drop=True)

    # published genotypes (Enbody et al. 2023), and GMM calls for everyone else
    pca["genotype"] = pca["IID"].map(published[locus])
    pca = classify_with_gmm(pca, GENOTYPES[locus])
    pca.to_csv(f"{OUTPUT}/genotypes_{locus}.csv", index=False)
    n_unknown = pca["prob_" + GENOTYPES[locus][0]].notna().sum()
    print(f"{locus}: {pca['genotype'].notna().sum()} published genotypes, {n_unknown} GMM calls "
          f"({(pca['confidence'] >= GENOTYPE_CONFIDENCE).sum() - pca['genotype'].notna().sum()} with confidence >= {GENOTYPE_CONFIDENCE})")

    # PCA colored by genotype and shaped by group: published genotypes plain, confident GMM calls outlined in black, uncertain ones in red
    # (predicted genotypes drawn over published ones, and within each, groups in GROUPS order so the founder's lineages and the founder are on top)
    fig, ax = plt.subplots(figsize=(9.5, 6))
    known, confident = pca["genotype"].notna(), pca["confidence"] >= GENOTYPE_CONFIDENCE
    for source, edge in [(known, "white"), (~known & confident, "black"), (~known & ~confident, "red")]:
        for group in GROUPS:
            subset = pca[source & (pca["Species"] == group)]
            if not subset.empty:
                sns.scatterplot(data=subset, x="PC1", y="PC2", hue="predicted_genotype", palette=GENOTYPE_COLORS, hue_order=GENOTYPES[locus],
                                marker=SPECIES_MARKERS[group], edgecolor=edge, alpha=0.7, legend=False, ax=ax)

    # legend: genotype colors, group shapes, and what the outlines mean
    blank = Line2D([0], [0], marker="", color="none")
    handles = [Line2D([0], [0], marker="", color="none", label="Genotype")]
    handles += [Line2D([0], [0], marker="o", color="w", markerfacecolor=GENOTYPE_COLORS[g], markersize=9, label=g) for g in GENOTYPES[locus]]
    handles += [blank, Line2D([0], [0], marker="", color="none", label="Group")]
    handles += [Line2D([0], [0], marker=m, color="w", markerfacecolor="lightgray", markeredgecolor="black", markersize=8, label=s) for s, m in SPECIES_MARKERS.items()]
    handles += [blank, Line2D([0], [0], marker="", color="none", label="Genotype source")]
    handles += [Line2D([0], [0], marker="o", color="w", markerfacecolor="lightgray", markeredgecolor=e, markersize=9, label=l)
                for e, l in [("white", "published"), ("black", f"predicted (≥ {GENOTYPE_CONFIDENCE})"), ("red", f"predicted (< {GENOTYPE_CONFIDENCE})")]]
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    for text in leg.get_texts():
        if text.get_text() in ["Genotype", "Group", "Genotype source"]:
            text.set_fontweight("bold")
        if text.get_text() in ["fortis", "scandens"]:
            text.set_fontstyle("italic")

    ax.set_xlabel(f"PC1 ({eigenval[0] / eigenval.sum() * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({eigenval[1] / eigenval.sum() * 100:.1f}%)")
    ax.set_title(f"Locus {locus}")
    sns.despine(top=True, right=True, ax=ax)
    plt.tight_layout()
    save_figure(f"pca_{locus}.png", fig)
