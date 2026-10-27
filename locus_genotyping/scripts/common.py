import importlib.util
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# paths are found from this file's location, so scripts work from any working directory
GENOTYPING = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(GENOTYPING)
DATA = f"{GENOTYPING}/data"
OUTPUT = f"{GENOTYPING}/output"
FIGURES = f"{GENOTYPING}/figures"

# reuse the bird metadata, IBD segments, and plot helpers from figures_ibis (loaded by path since all are named common.py)
_spec = importlib.util.spec_from_file_location("figures_common", f"{REPO}/figures_ibis/scripts/common.py")
sys.modules["figures_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["figures_common"])
from figures_common import (COLOR_MAP, FOCAL_SPECIES, FOUNDER, GROUPS, add_beak_body_ratio, add_bill_pcs, add_body_pc, expand_years_alive, format_year_axis,
                            founder_ibd_per_bird, load_chromosomes, load_master, load_trimmed_segments, set_plot_theme, style_legend)

# genotype calls from the GMM below this posterior probability are treated as missing
GENOTYPE_CONFIDENCE = 0.99

# a bird is in IBD at a locus if its trimmed segments cover at least this fraction of the locus (same 80% rule as the genomic windows)
LOCUS_IBD_COVERAGE = 0.8

# genotypes at each large-effect beak locus (Enbody et al. 2023): small / large alleles, or blunt / pointed at G07
GENOTYPES = {"G01": ["SS", "SL", "LL"], "G03": ["SS", "SL", "LL"], "G07": ["BB", "PB", "PP"],
             "G27": ["SS", "SL", "LL"], "G29": ["SS", "SL", "LL"], "G30": ["SS", "SL", "LL"]}


def load_loci():
    """The 6 large-effect beak loci (chrom, pos_start, pos_end) with short names (G01, ...), sorted by name."""
    loci = pd.read_csv(f"{DATA}/loci_positions.csv")
    loci["name"] = loci["locus"].str.split("_").str[-1]
    return loci.sort_values("name").reset_index(drop=True)


def locus_ibd_fraction(master):
    """Fraction of each locus covered by each bird's trimmed IBD segments, for every bird in the IBIS run (0 if none); birds x loci."""
    seg = load_trimmed_segments()
    seg["chrom_name"] = "chr" + seg["chrom"].map(load_chromosomes())
    birds = founder_ibd_per_bird(master)["ind"]
    fractions = {}
    for _, locus in load_loci().iterrows():
        segs = seg[seg["chrom_name"] == locus["chrom"]]
        overlap = (np.minimum(segs["phys_end_pos"] * 1e6, locus["pos_end"]) - np.maximum(segs["phys_start_pos"] * 1e6, locus["pos_start"])).clip(lower=0)
        fractions[locus["name"]] = birds.map(overlap.groupby(segs["ind"]).sum() / (locus["pos_end"] - locus["pos_start"])).fillna(0).to_numpy()
    return pd.DataFrame(fractions, index=birds.to_numpy())


def load_genotypes(locus):
    """Genotype of each bird at a locus: the published call if there is one, otherwise the GMM call if confident enough (else missing)."""
    calls = pd.read_csv(f"{OUTPUT}/genotypes_{locus}.csv").set_index("IID")
    return calls["predicted_genotype"].where(calls["confidence"] >= GENOTYPE_CONFIDENCE)


def save_figure(path, fig=None):
    """Save (and close) a figure with the same settings as figures_ibis (path relative to locus_genotyping/figures/)."""
    fig = fig or plt.gcf()
    fig.savefig(f"{FIGURES}/{path}", dpi=450, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
