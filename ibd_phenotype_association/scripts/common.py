import importlib.util
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# paths are found from this file's location, so scripts work from any working directory
ASSOCIATION = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(ASSOCIATION)
DATA = f"{ASSOCIATION}/data"
OUTPUT = f"{ASSOCIATION}/output"
RESULTS = f"{OUTPUT}/ibd_phenotype_association.csv"

# reuse the IBD, trait, window, and LMM helpers from figures_ibis (loaded by path since all are named common.py)
_spec = importlib.util.spec_from_file_location("figures_common", f"{REPO}/figures_ibis/scripts/common.py")
sys.modules["figures_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["figures_common"])
from figures_common import (COLOR_MAP, FOCAL_SPECIES, WINDOW_COVERAGE, add_beak_body_ratio, add_bill_pcs, add_body_pc, eigen_grm, lmm_lrt,
                            load_bim, load_chromosomes, load_founder_coef, load_loco_grm, load_master, load_trimmed_segments, make_windows,
                            set_plot_theme, window_coverage)

# traits tested, and the covariates of each (beak size / pointedness also control for body size, i.e. allometry)
TRAITS = ["PC(Body Size)", "PC(Beak Size)", "PC(Beak Pointedness)", "Beak:Body Ratio"]
COVARIATES = {"PC(Body Size)": ["sex", "loco_ibd"], "PC(Beak Size)": ["sex", "loco_ibd", "PC(Body Size)"],
              "PC(Beak Pointedness)": ["sex", "loco_ibd", "PC(Body Size)"], "Beak:Body Ratio": ["sex", "loco_ibd"]}
MEASUREMENTS = ["weight", "wing", "tarsus", "bill.length", "bill.depth", "bill.width"]

# a window is only tested if its IBD (test) group has at least this many birds, with at least MIN_PER_SEX of each sex
MIN_TEST = {"fortis": 15, "scandens": 10}
MIN_PER_SEX = 3
FDR = 0.05


def load_birds(spp):
    """Birds of a species in the IBIS run with all six measurements and known sex, with the body / beak size PCs and the beak:body ratio fit
    within the species."""
    master = load_master()
    in_ibis = master["meaningful.unique"].isin(load_founder_coef()["ind"])
    birds = master[(master["Species"] == spp) & in_ibis].dropna(subset=MEASUREMENTS + ["sex"]).rename(columns={"meaningful.unique": "ind"})
    return add_beak_body_ratio(add_bill_pcs(add_body_pc(birds.copy()))).reset_index(drop=True)


def loco_ibd_fractions():
    """Fraction of SNPs in (trimmed) IBD with the founder outside each chromosome, for every bird with IBD (bird x IBIS chromosome code)."""
    bim = load_bim()
    snps = bim.groupby("chrom").size()
    ibd = load_trimmed_segments().groupby(["ind", "chrom"])["SNP_count"].sum().unstack(fill_value=0).reindex(columns=snps.index, fill_value=0)
    return (ibd.sum(axis=1).to_numpy()[:, None] - ibd) / (len(bim) - snps)


def save_figure(path, fig=None):
    """Save (and close) a figure with the same settings as figures_ibis (path relative to ibd_phenotype_association/figures/)."""
    fig = fig or plt.gcf()
    fig.savefig(f"{ASSOCIATION}/figures/{path}", dpi=450, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
