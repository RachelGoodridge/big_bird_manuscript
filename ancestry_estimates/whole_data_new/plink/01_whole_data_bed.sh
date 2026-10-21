#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2

# conversion to plink1 .bed/.bim/.fam (for ADMIXTURE) plus a PCA for QC/visualization
$PLINK2 --pfile ../../training_new/plink/finches_pruned --allow-extra-chr --autosome-num 30 \
  --pca 10 --threads $THREADS --make-bed --out finches

echo "whole-data individuals:"
wc -l finches.fam
