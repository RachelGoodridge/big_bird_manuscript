import os
import pandas as pd
from common import CONTROLS, FOCAL_BIRDS, ORIGINAL, FOUNDER, load_master, read_bim
from ibis_common import VCF_SAMPLES, load_chr_map

# the Española conirostris in the vcf (the same birds as in the PCA and the ADMIXTURE panel of figure 1)
master = load_master(island=None)
with open(VCF_SAMPLES) as f:
    vcf_birds = {line.strip() for line in f}
espanola = master[(master["Island"] == "Espanola") & (master["Species"] == "conirostris") & master["meaningful.unique"].isin(vcf_birds)]
espanola = sorted(espanola["meaningful.unique"])
print(f"{len(espanola)} Española conirostris in the vcf: {', '.join(espanola)}")

# samples: every Daphne bird of the original run, and the Española conirostris; the founder and the Española conirostris come first (the
# order of focal_birds.txt), so that IBIS compares only them with every other bird (-setIndexStart / -setIndexEnd)
with open(f"{ORIGINAL}/dataset/Daphne_birds_for_ibis.txt") as f:
    daphne = [line.strip() for line in f if line.strip()]
assert FOUNDER in daphne and not set(espanola) & set(daphne)
focal = [FOUNDER] + espanola
samples = focal + [bird for bird in daphne if bird != FOUNDER]
os.makedirs(f"{CONTROLS}/dataset", exist_ok=True)
for path, birds in [(FOCAL_BIRDS, focal), (f"{CONTROLS}/dataset/samples.txt", samples)]:
    with open(path, "w") as f:
        f.writelines(f"{bird}\n" for bird in birds)
print(f"{len(samples)} samples ({len(daphne)} Daphne birds of the original run + {len(espanola)} Española conirostris)")

# SNPs: those of the original run's plink files (its downsampled SNPs that passed its missing data filter), by their vcf chromosome names
bim = read_bim(f"{ORIGINAL}/plink/dataset")
names = load_chr_map(f"{ORIGINAL}/dataset/chr_map.txt")  # integer code -> name without "chr"
snps = pd.DataFrame({"chrom": "chr" + bim["chrom"].map(names).str.replace("^chr", "", regex=True), "pos": bim["pos"]})
downsampled = pd.read_csv(f"{ORIGINAL}/dataset/downsampled_snps.tsv.gz", sep="\t", header=None, names=["chrom", "pos"])
assert snps.merge(downsampled).shape[0] == len(snps), "every SNP should be in the original downsampled SNP list"
snps.to_csv(f"{CONTROLS}/dataset/snps.tsv", sep="\t", header=False, index=False)
print(f"{len(snps)} SNPs of the original run")
