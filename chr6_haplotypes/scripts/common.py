import importlib.util
import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# paths are found from this file's location, so scripts work from any working directory
HAPLOTYPES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(HAPLOTYPES)
DATA = f"{HAPLOTYPES}/data"
WORK = f"{HAPLOTYPES}/work"  # plink files (large, not tracked)
OUTPUT = f"{HAPLOTYPES}/output"
FIGURES = f"{HAPLOTYPES}/figures"
CANDIDATES = f"{DATA}/candidates.tsv"  # birds passing the IBD and ancestry rules (select_birds.py)
BIRDS = f"{DATA}/birds.tsv"  # every bird analyzed and its group: the candidates passing the regional check (regional_check.py)
FST_GROUPS = f"{DATA}/fst_groups.tsv"  # the two FST groups (fst_groups.py)
FST_WINDOWS = f"{OUTPUT}/fst_windows.tsv"
PEAK_REGIONS = f"{OUTPUT}/fst_peak_regions.tsv"  # FST peaks (fst_vs_recombination.py)
VCF_SAMPLES = "/local/storage/all_birds_highcov_lowcov/samples.txt"  # birds in the phased VCF
GMAP = "/local/storage/gmap_files/chr6.gmap"  # chromosome 6 genetic map (columns pos, chr, cM)

# genome-wide ADMIXTURE ancestry (K = 9) of every bird in the VCF, and the columns of each species' clusters: fortis (V1), scandens (V6), conirostris
# (V2, the Española conirostris and the founder) and magnirostris (V3 + V8, the two magnirostris clusters, summed)
ADMIXTURE = f"{REPO}/ancestry_estimates/whole_data_new/admixture/finches"
ANCESTRY_CLUSTERS = {"fortis": [0], "scandens": [5], "conirostris": [1], "magnirostris": [2, 7]}
MIN_ANCESTRY = 0.5  # birds with less genome-wide ancestry of their own species are left out

# the chromosome 6 region and the IBD rule of the region test in ibd_phenotype_association (figure 4B), loaded by path since all are named common.py
_spec = importlib.util.spec_from_file_location("association_common", f"{REPO}/ibd_phenotype_association/scripts/common.py")
sys.modules["association_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["association_common"])
from association_common import CHR6_REGION, REGION_COVERAGE
from figures_common import load_chromosomes, load_founder_coef, load_master, load_trimmed_segments, set_plot_theme, window_coverage

CORE_BIRDS = 100  # reference birds of each species for the regional check (below)

# groups of birds compared (in plotting order): label, color (the project's species colors), island, species (as in load_master), the species whose
# ancestry is checked against MIN_ANCESTRY, and the most birds kept, those with the most ancestry (None: all). The magnirostris are the 39 with the most
# magnirostris ancestry, so that the donors (below) number 50
GROUPS = pd.DataFrame([
    ("fortis",           "fortis (no IBD)",                         "#1C86EE", "Daphne",   "fortis",       "fortis",       CORE_BIRDS),
    ("scandens_no_ibd",  "scandens (no IBD)",                       "#43CD80", "Daphne",   "scandens",     "scandens",     None),
    ("scandens_ibd",     f"scandens (IBD ≥ {REGION_COVERAGE:.0%})", "#1B6B3A", "Daphne",   "scandens",     "scandens",     None),
    ("conirostris_esp",  "conirostris (Española)",                  "#FF3030", "Espanola", "conirostris",  "conirostris",  None),
    ("founder",          "founder",                                 "#FF3030", "Daphne",   "founder",      "conirostris",  None),
    ("magnirostris_dap", "magnirostris (Daphne)",                   "#EEB422", "Daphne",   "magnirostris", "magnirostris", 39),
], columns=["group", "label", "color", "island", "species", "ancestry", "max_birds"]).set_index("group")

# the two groups compared (FST scan and diagnostic SNPs): conirostris, the founder, and magnirostris (the donors) vs scandens without IBD. Alleles are
# counted as copies of the donor allele (the allele more common in the donors); SNPs whose frequency differs by at least DIAGNOSTIC between the two
# groups tell the two kinds of haplotype apart
DONORS = ["conirostris_esp", "founder", "magnirostris_dap"]
DIAGNOSTIC = 0.8

# regional check: scandens without IBD left out when, in the chromosome 6 region, they are closer (mean genotype distance, LD-pruned SNPs) to the
# reference fortis than to the reference scandens; the reference birds are the CORE_BIRDS of each with the most ancestry. The fortis are used only
# as this reference (REFERENCE_ONLY), not analyzed
REGIONAL_CHECK = {"scandens_no_ibd": "fortis"}
REFERENCE_ONLY = ["fortis"]

CHROM, START, END = CHR6_REGION


def load_birds():
    return pd.read_csv(BIRDS, sep="\t")


def region_ibd(birds):
    """Fraction of the chromosome 6 region covered by each bird's (summed, trimmed) IBD segments with the founder (0 if none)."""
    chrom_code = {name: code for code, name in load_chromosomes().items()}[CHROM]
    seg = load_trimmed_segments()
    region = pd.DataFrame({"chrom": [chrom_code], "start": [START], "end": [END]})
    coverage = window_coverage(seg[seg["ind"].isin(birds)], region).set_index("ind")["coverage"]
    return coverage.reindex(birds).fillna(0).to_numpy()


def species_ancestry(species):
    """Genome-wide ancestry from a species' clusters (ADMIXTURE, K = 9, summed over ANCESTRY_CLUSTERS) of every bird in the VCF."""
    q = pd.read_csv(f"{ADMIXTURE}.9.Q", sep=r"\s+", header=None)
    ids = pd.read_csv(f"{ADMIXTURE}.fam", sep=r"\s+", header=None)[1]
    return pd.Series(q[ANCESTRY_CLUSTERS[species]].sum(axis=1).to_numpy(), index=ids)


def load_genotypes():
    """Genotypes of every SNP in the region (work/region plink files): SNP positions, bird IDs, and a birds x SNPs array of copies of one of the
    SNP's alleles (0 / 1 / 2, NaN if missing; which allele does not matter, as donor_alleles finds the donor one)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)  # pandas_plink imports pandera the deprecated way
        from pandas_plink import read_plink
        bim, fam, bed = read_plink(f"{WORK}/region", verbose=False)
    return bim["pos"].to_numpy(), fam["iid"].to_numpy(), bed.compute().T


def donor_alleles(genotypes, groups):
    """For each SNP (genotypes: birds x SNPs of allele counts, groups: each bird's group), whether the counted allele is the donor allele, and whether
    the SNP is diagnostic (frequency differing by >= DIAGNOSTIC between the donors and the scandens without IBD, the two FST groups)."""
    freq_donor = np.nanmean(genotypes[np.isin(groups, DONORS)], axis=0) / 2
    freq_sca = np.nanmean(genotypes[groups == "scandens_no_ibd"], axis=0) / 2
    return freq_donor >= freq_sca, np.abs(freq_donor - freq_sca) >= DIAGNOSTIC


def save_figure(path, fig=None):
    """Save (and close) a figure as a png with the same settings as figures_ibis (path relative to chr6_haplotypes/figures/)."""
    fig = fig or plt.gcf()
    fig.savefig(f"{FIGURES}/{path}", dpi=450, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
