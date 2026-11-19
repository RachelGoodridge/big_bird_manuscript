import numpy as np
import pandas as pd
from common import FOUNDER, IBIS_PREFIX, ORIGINAL, SEG_COLUMNS, load_focal_birds

# the founder's segments and kinship coefficients with the Daphne birds should be identical in both runs, as IBIS compares each pair on its
# own and both runs have the same SNPs, genetic map and parameters
seg = pd.read_csv(f"{IBIS_PREFIX}.seg", sep="\t", header=None, names=SEG_COLUMNS)
seg = seg[(seg["sample1"] == FOUNDER) | (seg["sample2"] == FOUNDER)].copy()
seg["ind"] = np.where(seg["sample1"] == FOUNDER, seg["sample2"], seg["sample1"])
coef = pd.read_csv(f"{IBIS_PREFIX}.coef", sep="\t")
coef = coef[(coef["Individual1"] == FOUNDER) | (coef["Individual2"] == FOUNDER)].copy()
coef["ind"] = np.where(coef["Individual1"] == FOUNDER, coef["Individual2"], coef["Individual1"])

original_seg = pd.read_pickle(f"{ORIGINAL}/output/seg_all_founder.pkl")
original_coef = pd.read_pickle(f"{ORIGINAL}/output/coef_all_founder.pkl")
espanola = load_focal_birds()[1:]

cols = ["ind", "chrom", "phys_start_pos", "phys_end_pos", "IBD_type", "genetic_start_pos", "genetic_end_pos", "marker_count", "error_count"]
sort = lambda df: df[cols].sort_values(cols).reset_index(drop=True)
new = sort(seg[~seg["ind"].isin(espanola)])
assert new.equals(sort(original_seg)), "the founder's segments differ from the original run"
print(f"the founder's {len(new)} segments with Daphne birds are identical to the original run's")

coef_cols = [c for c in original_coef.columns if c not in ("Individual1", "Individual2")]
sort = lambda df: df[coef_cols].sort_values("ind").reset_index(drop=True)
new_coef = sort(coef[~coef["ind"].isin(espanola)])
assert new_coef.equals(sort(original_coef)), "the founder's kinship coefficients differ from the original run"
print(f"the founder's {len(new_coef)} kinship coefficients with Daphne birds are identical to the original run's")

# the founder's relatedness to the Española conirostris (not in the original run)
kin = coef[coef["ind"].isin(espanola)].set_index("ind")[["Kinship_Coefficient", "Degree"]].reindex(espanola)
print("founder vs. each Española conirostris:\n" + kin.to_string())
