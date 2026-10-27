import pandas as pd
from common import DATA, load_loci

# SNPs genotyped by Enbody et al. (2023) that fall within each large-effect beak locus, one bcftools regions file (CHROM, POS) per locus
snps = pd.read_csv(f"{DATA}/autosomes_Daphne_cluster1_multivariate_PCsp_lmm4.bin.pos.genos", sep="\t", usecols=["CHROM", "POS"])
for _, locus in load_loci().iterrows():
    inside = snps[(snps["CHROM"] == locus["chrom"]) & snps["POS"].between(locus["pos_start"], locus["pos_end"])]
    inside.to_csv(f"work/regions_{locus['name']}.txt", sep="\t", index=False, header=False)
    print(f"{locus['name']}: {len(inside)} SNPs")
