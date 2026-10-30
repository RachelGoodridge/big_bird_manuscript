import importlib.util
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

# paths are found from this file's location, so scripts work from any working directory
DROUGHT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(DROUGHT)
DATA = f"{DROUGHT}/data"

# reuse the IBD, bird metadata and plot helpers from figures_ibis (loaded by path since all are named common.py)
_spec = importlib.util.spec_from_file_location("figures_common", f"{REPO}/figures_ibis/scripts/common.py")
sys.modules["figures_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["figures_common"])
from figures_common import (COLOR_MAP, GENOTYPING, IBD_LABEL, founder_ibd_per_bird, load_chromosomes, load_master, load_trimmed_segments,
                            set_plot_theme, style_legend)


# load the drought survival list (Erik's published table) for Daphne fortis, with IBD and beak size
def load_drought_fortis(master):
    drought_surv = pd.read_csv(f"{DATA}/table_of_drought_fortis.csv").set_index("meaningful.unique")["drought_surv"]
    birds = master[master["meaningful.unique"].isin(drought_surv.index)].copy()
    birds["outcome"] = birds["meaningful.unique"].map(drought_surv)
    birds = birds[birds["Species"] == "fortis"].copy()
    ibd = founder_ibd_per_bird(master).set_index("ind")["fraction_snps_ibd"]
    birds["fraction_snps_ibd"] = birds["meaningful.unique"].map(ibd).fillna(0)
    birds["has_ibd"] = birds["fraction_snps_ibd"] > 0
    birds["beak_sum"] = birds["bill.length"] + birds["bill.depth"] + birds["bill.width"]
    return birds


def save_figure(path, fig=None):
    """Save (and close) a figure with the same settings as figures_ibis (path relative to drought_survivorship/figures/)."""
    fig = fig or plt.gcf()
    fig.savefig(f"{DROUGHT}/figures/{path}", dpi=450, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
