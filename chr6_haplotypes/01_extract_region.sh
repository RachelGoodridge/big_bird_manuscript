#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data work output figures

VCF=/local/storage/all_birds_highcov_lowcov/autosomes_2022chrZ_all_highcov_lowcov_phased.vcf.gz
PLINK2=/programs/plink2_linux_20230721/plink2
REGION=chr6:10000000-30500000  # CHR6_REGION in ibd_phenotype_association/scripts/common.py

# candidate birds of each group (IBD and ancestry rules, data/candidates.tsv)
python scripts/select_birds.py

# the chromosome 6 region for them, keeping SNPs genotyped in >= 95% of them (as in tune_trim_ibis) with MAF >= 1%, and its LD-pruned SNPs
# (r2 < 0.2 in 50 SNP windows, step 10)
bcftools view --regions $REGION --samples-file data/candidates.txt --threads 8 -Oz -o work/candidates.vcf.gz $VCF
$PLINK2 --vcf work/candidates.vcf.gz --allow-extra-chr --autosome-num 30 --set-missing-var-ids @:# --geno 0.05 --maf 0.01 \
    --make-bed --out work/candidates --silent
rm work/candidates.vcf.gz
$PLINK2 --bfile work/candidates --indep-pairwise 50 10 0.2 --out work/candidates_pruned --silent
$PLINK2 --bfile work/candidates --extract work/candidates_pruned.prune.in --make-bed --out work/candidates_pruned --silent

# the birds analyzed: fortis / scandens closer to the other species in the region left out (data/birds.tsv), and the region for them
python scripts/regional_check.py
$PLINK2 --bfile work/candidates --keep data/birds.txt --geno 0.05 --maf 0.01 --make-bed --out work/region --silent
