#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2

# drop variants missing in > 10% of samples
$PLINK2 --pfile finches_00_raw --allow-extra-chr --autosome-num 30 \
  --geno 0.1 --threads $THREADS --make-pgen --out finches_01_geno
