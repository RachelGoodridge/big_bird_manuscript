import shutil
import pandas as pd
from common import CONTROLS, ORIGINAL, load_focal_birds, read_bim
from ibis_common import load_chr_map

# the new plink files have the original run's SNPs, in the same order and with the same alleles (plink names chromosomes from the vcf; the
# original run replaced them with integer codes, see tune_trim_ibis/scripts/make_map.py)
new = pd.read_csv(f"{CONTROLS}/plink/dataset.bim", sep="\t", header=None, names=["chrom", "snp", "cm", "pos", "a0", "a1"], dtype={"chrom": str})
original = read_bim(f"{ORIGINAL}/plink/dataset")
names = load_chr_map(f"{ORIGINAL}/dataset/chr_map.txt")
codes = {name.replace("chr", ""): code for code, name in names.items()}
assert len(new) == len(original), f"{len(new)} SNPs, the original run has {len(original)}"
assert (new["chrom"].str.replace("chr", "").map(codes).to_numpy() == original["chrom"].to_numpy()).all()
for col in ["pos", "a0", "a1"]:
    assert (new[col].to_numpy() == original[col].to_numpy()).all(), f"{col} differs from the original run"

# so the original bim file (integer chromosome codes, SNP names and genetic map) describes the new plink files exactly
shutil.copyfile(f"{ORIGINAL}/plink/dataset.bim", f"{CONTROLS}/plink/dataset.bim")
print(f"the {len(new)} SNPs match the original run (positions and alleles); copied its bim file (chromosome codes and genetic map)")

# the founder and the Española conirostris come first, as IBIS needs
fam = pd.read_csv(f"{CONTROLS}/plink/dataset.fam", sep=r"\s+", header=None)
focal = load_focal_birds()
assert list(fam[1][:len(focal)]) == focal, "the founder and the Española conirostris should be the first samples"
print(f"{len(fam)} samples, the first {len(focal)} the founder and the Española conirostris")

# missing data of the Española birds (the SNPs were filtered for missing data in the Daphne birds only)
miss = pd.read_csv(f"{CONTROLS}/plink/dataset.smiss", sep="\t").set_index("IID")["F_MISS"]
print("fraction of SNPs missing: founder %.4f, Española conirostris %s, Daphne median %.4f"
      % (miss[focal[0]], ", ".join(f"{miss[b]:.4f}" for b in focal[1:]), miss.drop(focal).median()))
