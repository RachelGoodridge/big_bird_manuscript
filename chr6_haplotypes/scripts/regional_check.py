import warnings
import numpy as np
import pandas as pd
from common import BIRDS, CANDIDATES, CORE_BIRDS, DATA, GROUPS, REFERENCE_ONLY, REGIONAL_CHECK, WORK

# the birds analyzed: the candidates (select_birds.py) minus the reference-only fortis and the scandens without IBD that, in the chromosome 6 region,
# are closer to the fortis than to the other scandens. Distance between two birds: mean |difference in allele counts| / 2 over the LD-pruned SNPs (work/candidates_pruned)
# both are genotyped at, i.e. 1 - IBS. Each bird's distance to a species is its mean distance to that species' reference birds (the CORE_BIRDS with the
# most genome-wide ancestry, not counting itself)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)  # pandas_plink imports pandera the deprecated way
    from pandas_plink import read_plink
    _, fam, bed = read_plink(f"{WORK}/candidates_pruned", verbose=False)
genotypes = bed.compute().T.astype(np.float32)  # birds x SNPs
candidates = pd.read_csv(CANDIDATES, sep="\t").set_index("ind").loc[fam["iid"]].reset_index()


def mean_distance(rows, ref):
    """Mean distance of each bird in rows to the birds in ref (a bird's distance to itself not counted)."""
    total, count = np.zeros(len(rows)), np.zeros(len(rows))
    for j in ref:
        d = np.nanmean(np.abs(genotypes[rows] - genotypes[j]), axis=1) / 2
        d[rows == j] = np.nan
        total += np.nan_to_num(d)
        count += ~np.isnan(d)
    return total / count


core = {g: candidates[candidates["group"] == g].nlargest(CORE_BIRDS, "ancestry").index.to_numpy() for g in {*REGIONAL_CHECK, *REGIONAL_CHECK.values()}}
checks = []
for g, other in REGIONAL_CHECK.items():
    rows = candidates.index[candidates["group"] == g].to_numpy()
    check = pd.DataFrame({"ind": candidates.loc[rows, "ind"], "group": g, "ancestry": candidates.loc[rows, "ancestry"],
                          "dist_own": mean_distance(rows, core[g]), "dist_other": mean_distance(rows, core[other])})
    check["closer_to_other"] = check["dist_other"] < check["dist_own"]
    print(f"{g}: {check['closer_to_other'].sum()} of {len(check)} closer to the {GROUPS.loc[other, 'species']} than to the {GROUPS.loc[g, 'species']} "
          f"in the region (left out)")
    checks.append(check)
checks = pd.concat(checks)
checks.sort_values(["group", "closer_to_other", "ancestry"]).to_csv(f"{DATA}/regional_check.tsv", sep="\t", index=False, float_format="%.4f")

birds = candidates[~candidates["ind"].isin(checks.loc[checks["closer_to_other"], "ind"]) & ~candidates["group"].isin(REFERENCE_ONLY)]
birds.to_csv(BIRDS, sep="\t", index=False)
birds["ind"].to_csv(f"{DATA}/birds.txt", index=False, header=False)
print(birds["group"].value_counts().reindex(GROUPS.index.drop(REFERENCE_ONLY)).to_string())
print(f"saved {len(birds)} birds to {BIRDS}")
