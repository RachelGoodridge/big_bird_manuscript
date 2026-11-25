import importlib.util
import os
import sys

# Española conirostris as controls for the founder: an IBIS run comparing the founder and each Española conirostris with every Daphne bird of
# the original run (tune_trim_ibis), made as comparable to it as possible: the same SNPs (those passing its missing data filter), the same genetic
# map, the same IBIS parameters, the same IBD threshold (its negative control max error) and the same trimmed false positive regions. Nothing in
# tune_trim_ibis is written; its files are only read. Paths are found from this file's location, so scripts work from any working directory
CONTROLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(CONTROLS)
ORIGINAL = f"{REPO}/tune_trim_ibis"  # the original IBIS run

FULL_VCF = "/local/storage/all_birds_highcov_lowcov/autosomes_chrZ_all_highcov_lowcov.vcf.gz"
IBIS_PREFIX = f"{CONTROLS}/output/espanola"  # IBIS output (.seg, .coef)
FOCAL_BIRDS = f"{CONTROLS}/dataset/focal_birds.txt"  # the founder and the Española conirostris, first in the plink files

# the original run's helpers (loaded by path since the files are all named common.py)
if "ibis_common" not in sys.modules:
    _spec = importlib.util.spec_from_file_location("ibis_common", f"{ORIGINAL}/scripts/common.py")
    sys.modules["ibis_common"] = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(sys.modules["ibis_common"])
from ibis_common import FOUNDER, SEG_COLUMNS, count_snps, load_master, negative_controls, read_bed, read_bim, subtract_regions  # noqa: E402


def load_focal_birds():
    """The founder, then the Española conirostris (the order of the plink files)."""
    with open(FOCAL_BIRDS) as f:
        return [line.strip() for line in f if line.strip()]


def original_threshold():
    """The original run's IBD threshold: the largest fraction of SNPs any negative control shares in IBD with the founder (computed as in
    tune_trim_ibis/scripts/control_segments.py)."""
    import pandas as pd
    coef = pd.read_pickle(f"{ORIGINAL}/output/coef_all_founder.pkl")
    neg_ids = set(negative_controls(load_master())["meaningful.unique"]) & set(coef["ind"])
    n_snps = len(read_bim(f"{ORIGINAL}/plink/dataset"))
    founder_frac = pd.read_pickle(f"{ORIGINAL}/output/seg_all_founder.pkl").groupby("ind")["SNP_count"].sum() / n_snps
    return founder_frac.reindex(sorted(neg_ids)).fillna(0).max()


############################ minimum segment length (a test, see 05_min_length.sh) ############################
# segments must be longer than the longest one the founder shares with any negative control, after each pair's segments separated by gaps of at
# most MERGE_GAP_CM are merged (IBIS often splits one long IBD segment into shorter pieces)
MERGE_GAP_CM = 1
MIN_LENGTH = f"{CONTROLS}/min_length_merge{MERGE_GAP_CM:g}cM"  # its segments and plots


def merged_length(seg, gap_cm=MERGE_GAP_CM):
    """Genetic length (cM) of the block each segment belongs to: each pair's consecutive segments on a chromosome are merged when the gap between
    them is at most gap_cm cM, from the first piece's start to the last one's end. seg needs columns focal, partner, chrom, genetic_start_pos and
    genetic_end_pos."""
    s = seg.sort_values(["focal", "partner", "chrom", "genetic_start_pos"])
    same = ((s["focal"] == s["focal"].shift()) & (s["partner"] == s["partner"].shift()) & (s["chrom"] == s["chrom"].shift())).to_numpy()
    close = (s["genetic_start_pos"] - s["genetic_end_pos"].shift() <= gap_cm).to_numpy()
    block = (~(same & close)).cumsum()
    start = s.groupby(block)["genetic_start_pos"].transform("min")
    end = s.groupby(block)["genetic_end_pos"].transform("max")
    return (end - start).reindex(seg.index)


def original_min_length():
    """The minimum segment length tested (cM): the longest merged block (see merged_length) the founder shares with any negative control in the
    original run (the segment counterpart of its threshold, the negative control max error). Segments must be in longer blocks."""
    import pandas as pd
    seg = pd.read_pickle(f"{ORIGINAL}/output/seg_all_founder.pkl").rename(columns={"ind": "partner"}).assign(focal=FOUNDER)
    return merged_length(seg)[seg["partner"].isin(negative_controls(load_master())["meaningful.unique"])].max()


def focal_pairs(seg, focal_birds):
    """Segments between each focal bird and every other bird, as columns focal and partner (pairs of two focal birds are kept once for each)."""
    import pandas as pd
    return pd.concat([seg[seg[focal].isin(focal_birds)].rename(columns={focal: "focal", partner: "partner"})
                      for focal, partner in [("sample1", "sample2"), ("sample2", "sample1")]], ignore_index=True)


def filter_and_trim(pairs, bim, threshold, exclude, min_cm=None):
    """The original run's filtering of focal pairs (tune_trim_ibis/scripts/control_segments.py): keep pairs sharing a larger fraction of SNPs in
    IBD than threshold, then (if min_cm is given) only their segments in merged blocks longer than min_cm cM (see merged_length), then trim the
    excluded regions, recount SNPs, drop pieces with 1 or fewer SNPs and convert positions to Mb (the trimmed pieces are returned)."""
    pairs = pairs.copy()
    pairs["SNP_count"] = count_snps(pairs, bim)
    fraction_snps_ibd = pairs.groupby(["focal", "partner"], observed=True)["SNP_count"].transform("sum") / len(bim)
    kept = pairs[fraction_snps_ibd > threshold]
    if min_cm is not None:
        kept = kept[merged_length(kept) > min_cm]
    kept = kept.reset_index(drop=True)
    trim_seg = subtract_regions(kept[["focal", "partner", "chrom", "phys_start_pos", "phys_end_pos"]], exclude)
    trim_seg["SNP_count"] = count_snps(trim_seg, bim)
    trim_seg = trim_seg[trim_seg["SNP_count"] > 1].copy()
    trim_seg["phys_start_pos"] /= 1e6
    trim_seg["phys_end_pos"] /= 1e6
    return trim_seg.reset_index(drop=True)
