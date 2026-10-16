#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2

# drop variants with minor allele frequency < 5% across the combined sample set;
# adjust --maf if you want stricter/looser. Variants are already biallelic so no
# --max-alleles filter is needed.
$PLINK2 --pfile finches_01_geno --allow-extra-chr --autosome-num 30 \
  --maf 0.05 --threads $THREADS --make-pgen --out finches_02_maf
