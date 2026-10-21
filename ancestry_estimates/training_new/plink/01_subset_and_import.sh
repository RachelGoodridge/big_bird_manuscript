#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2
FULL_VCF=/local/storage/all_birds_highcov_lowcov/autosomes_2022chrZ_all_highcov_lowcov_phased.vcf.gz

# chrZ is excluded: mostly missing and not used in IBIS
$PLINK2 --vcf $FULL_VCF --not-chr chrZ --allow-extra-chr --autosome-num 30 \
  --set-missing-var-ids @:# --threads $THREADS --make-pgen --out finches_00_raw
