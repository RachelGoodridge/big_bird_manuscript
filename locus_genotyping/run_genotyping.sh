#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p work figures

FULL_VCF=/local/storage/all_birds_highcov_lowcov/autosomes_2022chrZ_all_highcov_lowcov_phased.vcf.gz
PLINK2=/programs/plink2_linux_20230721/plink2

# the Enbody et al. (2023) SNPs within each large-effect beak locus
python scripts/locus_snps.py

# genotypes of every bird in the full vcf at those SNPs, and a PCA of them for each locus
for REGIONS in work/regions_*.txt; do
  LOCUS=$(basename $REGIONS .txt | sed 's/regions_//')
  bcftools view -R $REGIONS $FULL_VCF -Oz -o work/$LOCUS.vcf.gz
  bcftools index -f work/$LOCUS.vcf.gz
  $PLINK2 --vcf work/$LOCUS.vcf.gz --allow-extra-chr --set-missing-var-ids @:# --pca 10 --out work/pca_$LOCUS
done

# genotype calls (output/genotypes_G*.csv) and PCA plots
python scripts/call_genotypes.py

# figures
python scripts/plot_allele_freqs.py
python scripts/plot_genotype_vs_locus_ibd.py
