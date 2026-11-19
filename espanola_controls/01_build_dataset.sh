#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p dataset plink

FULL_VCF=/local/storage/all_birds_highcov_lowcov/autosomes_2022chrZ_all_highcov_lowcov_phased.vcf.gz
VCF=dataset/espanola_downsampled
PLINK2=/programs/plink2_linux_20230721/plink2

# the samples (the founder and the Española conirostris first, then every other Daphne bird of the original run) and the SNPs of the
# original run (tune_trim_ibis/plink/dataset.bim)
python scripts/select_samples.py

# get those birds and SNPs from the full vcf (as in tune_trim_ibis/01_build_dataset.sh)
bcftools view --threads 16 -S dataset/samples.txt -T dataset/snps.tsv $FULL_VCF -Oz -o $VCF.vcf.gz
bcftools index $VCF.vcf.gz
bcftools stats $VCF.vcf.gz > $VCF.stats

# convert to plink in the order of samples.txt. No missing data filter: the SNPs are already those that passed the original run's filter
# (refiltering with the Española birds added could drop a few), so both runs use exactly the same SNPs
$PLINK2 --vcf $VCF.vcf.gz --allow-extra-chr --autosome-num 95 --indiv-sort file dataset/samples.txt --make-bed --out plink/dataset
$PLINK2 --bfile plink/dataset --allow-extra-chr --autosome-num 95 --missing sample-only --out plink/dataset

# check that the SNPs match the original run's, and give them its chromosome codes and genetic map
python scripts/make_map.py

# clear up the large vcf file
rm $VCF.vcf.gz $VCF.vcf.gz.csi
