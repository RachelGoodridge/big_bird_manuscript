#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p dataset plink

FULL_VCF=/local/storage/all_birds_highcov_lowcov/autosomes_2022chrZ_all_highcov_lowcov_phased.vcf.gz
SAMPLES=dataset/Daphne_birds_for_ibis.txt
SNPS=dataset/downsampled_snps.tsv.gz
NEW_VCF=dataset/Daphne_birds_downsampled
PLINK2=/programs/plink2_linux_20230721/plink2

# list all Daphne birds in the vcf and the birds for tuning
python scripts/select_samples.py

# get all the birds from the large dataset and downsample to the SNPs in the genetic map
bcftools view --threads 16 -S $SAMPLES -T $SNPS $FULL_VCF -Oz -o $NEW_VCF.vcf.gz
bcftools index $NEW_VCF.vcf.gz
bcftools stats $NEW_VCF.vcf.gz > $NEW_VCF.stats

# then use plink to filter out SNPs with >5% missing data
$PLINK2 --vcf $NEW_VCF.vcf.gz --geno 0.05 --allow-extra-chr --autosome-num 95 --make-bed --out plink/dataset

# next add the genetic map to the bim file
python scripts/make_map.py

# clear up the massive vcf file
rm $NEW_VCF.vcf.gz $NEW_VCF.vcf.gz.csi