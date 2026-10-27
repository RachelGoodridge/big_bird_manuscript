import importlib.util
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import seaborn as sns

# paths are found from this file's location, so scripts work from any working directory
EPISTASIS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(EPISTASIS)

# reuse the genotype calls and locus IBD from locus_genotyping (loaded by path since all are named common.py)
_spec = importlib.util.spec_from_file_location("genotyping_common", f"{REPO}/locus_genotyping/scripts/common.py")
sys.modules["genotyping_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["genotyping_common"])
from genotyping_common import (GENOTYPES, LOCUS_IBD_COVERAGE, add_beak_body_ratio, add_bill_pcs, add_body_pc, load_genotypes, load_loci, load_master, locus_ibd_fraction,
                               set_plot_theme)

# traits compared between genotype combinations, with short names for file names
TRAITS = {"PC(Body Size)": "body_size", "PC(Beak Size)": "beak_size", "Beak:Body Ratio": "beak_body_ratio"}
SPECIES_NAMES = {"fortis", "scandens"}


def load_birds():
    """Daphne birds with body / beak size PCs and the beak:body ratio, and their genotype at each locus (missing if not confidently called)."""
    # traits fit once across all Daphne birds with every measurement (except magnirostris), so every group shares the same axes
    master = load_master()
    birds = master[master["Species"] != "magnirostris"].dropna(subset=["weight", "wing", "tarsus", "bill.length", "bill.depth", "bill.width"]).copy()
    birds = add_beak_body_ratio(add_bill_pcs(add_body_pc(birds)))

    birds = birds.set_index("meaningful.unique")
    for locus in GENOTYPES:
        birds[locus] = load_genotypes(locus)
    return birds


def split_by_ibd(birds, loci):
    """Add "-IBD" to the genotype of birds in IBD with the founder at the locus (birds outside the IBIS run count as not in IBD)."""
    in_ibd = locus_ibd_fraction(load_master()) >= LOCUS_IBD_COVERAGE
    birds = birds.copy()
    for locus in loci:
        ibd = birds.index.map(in_ibd[locus]).fillna(False).astype(bool)
        birds[locus] = birds[locus].mask(birds[locus].notna() & ibd, birds[locus] + "-IBD")
    return birds


def genotype_order(locus, by_ibd):
    return GENOTYPES[locus] + [f"{g}-IBD" for g in GENOTYPES[locus]] if by_ibd else GENOTYPES[locus]


def cell_summary(birds, locus1, locus2, trait, by_ibd):
    """Mean, standard deviation, and number of birds for each genotype combination (locus1 rows x locus2 columns)."""
    grouped = birds.dropna(subset=[locus1, locus2, trait]).groupby([locus1, locus2])[trait]
    tables = [grouped.agg(stat).unstack(locus2) for stat in ["mean", "std", "count"]]
    rows = [g for g in genotype_order(locus1, by_ibd) if g in tables[0].index]
    cols = [g for g in genotype_order(locus2, by_ibd) if g in tables[0].columns]
    return [t.reindex(index=rows, columns=cols) for t in tables]


def group_title(group):
    return rf"$\it{{{group}}}$" if group in SPECIES_NAMES else group


def plot_mean_heatmap(means, ns, locus1, locus2, trait, title, path):
    """Mean trait value (and number of birds) of each genotype combination, with the homozygous small genotypes at the bottom left."""
    fig, ax = plt.subplots(figsize=(max(4, len(means.columns) * 1.5 + 1.5), max(3.5, len(means.index) * 1.0 + 1.5)))
    annot = means.map(lambda v: "—" if pd.isna(v) else f"{v:.2f}") + ns.map(lambda n: "" if pd.isna(n) else f"\n(n={int(n)})")
    values = means.astype(float).to_numpy()
    norm = mcolors.TwoSlopeNorm(vmin=np.nanmin(values), vcenter=np.nanmean(values), vmax=np.nanmax(values))
    sns.heatmap(means.astype(float), annot=annot, fmt="", cmap="RdYlBu_r", norm=norm, linewidths=0.5, linecolor="white",
                cbar_kws={"label": trait}, ax=ax)
    ax.grid(False)
    ax.invert_yaxis()
    ax.tick_params(axis="y", rotation=0)
    ax.set_xlabel(locus2)
    ax.set_ylabel(locus1)
    ax.set_title(title, pad=10)
    plt.tight_layout()
    save_figure(path, fig)


def save_figure(path, fig=None):
    """Save (and close) a figure with the same settings as figures_ibis (path relative to epistasis/figures/)."""
    fig = fig or plt.gcf()
    fig.savefig(f"{EPISTASIS}/figures/{path}", dpi=450, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
