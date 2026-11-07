import numpy as np
import pandas as pd
from common import END, FST_WINDOWS, START, WORK

WINDOW = 100_000
STEP = 50_000
MIN_SNPS = 50  # windows with fewer SNPs are skipped


########################### PER-SNP FST ###########################
fst = pd.read_csv(f"{WORK}/fst.con_mag.scandens.fst.var", sep="\t")
fst.columns = fst.columns.str.lstrip("#")

# Hudson FST numerator and denominator per SNP from each group's allele counts (the estimator plink2 uses, Bhatia et al. 2013), so that windows
# are a ratio of averages
counts = {g: pd.read_csv(f"{WORK}/freq_{g}.acount", sep="\t").set_index("ID") for g in ["con_mag", "scandens"]}
n1, n2 = (fst["ID"].map(counts[g]["OBS_CT"]) for g in counts)
p1, p2 = (fst["ID"].map(counts[g]["ALT_CTS"]) / n for g, n in zip(counts, [n1, n2]))
fst["NUM"] = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
fst["DEN"] = p1 * (1 - p2) + p2 * (1 - p1)
fst = fst.dropna(subset=["HUDSON_FST", "NUM", "DEN"]).sort_values("POS").reset_index(drop=True)
print(f"{len(fst)} SNPs with a defined FST; FST over the region (ratio of averages): {fst['NUM'].sum() / fst['DEN'].sum():.3f}")

########################### WINDOWS ###########################
pos = fst["POS"].to_numpy()
cum_num = np.concatenate([[0], np.cumsum(fst["NUM"])])
cum_den = np.concatenate([[0], np.cumsum(fst["DEN"])])
windows = []
for start in range(START, END - WINDOW + 1, STEP):
    lo, hi = np.searchsorted(pos, start), np.searchsorted(pos, start + WINDOW)
    if hi - lo >= MIN_SNPS:
        windows.append({"start": start, "end": start + WINDOW, "n_snps": hi - lo, "fst": (cum_num[hi] - cum_num[lo]) / (cum_den[hi] - cum_den[lo])})
windows = pd.DataFrame(windows)
windows["mid"] = (windows["start"] + windows["end"]) / 2
windows["z_fst"] = (windows["fst"] - windows["fst"].mean()) / windows["fst"].std()
windows.to_csv(FST_WINDOWS, sep="\t", index=False, float_format="%.5f")
print(f"{len(windows)} windows with >= {MIN_SNPS} SNPs, highest FST {windows['fst'].max():.3f} (Z = {windows['z_fst'].max():.2f})")
