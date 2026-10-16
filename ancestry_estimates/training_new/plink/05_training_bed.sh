#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2

# 1) build a KING kinship matrix among the 156 candidate training individuals
$PLINK2 --pfile finches_pruned --keep keep_list.txt --allow-extra-chr --autosome-num 30 \
  --make-king triangle bin --threads $THREADS --out finches

# 2) drop one member of every related pair above the 0.177 (~3rd-degree) cutoff,
#    then PCA (for QC/visualization) and export to plink1 .bed/.bim/.fam for ADMIXTURE
$PLINK2 --pfile finches_pruned --keep keep_list.txt --allow-extra-chr --autosome-num 30 \
  --king-cutoff finches 0.177 --pca 10 --threads $THREADS --make-bed --out finches

echo "training individuals remaining after relatedness filtering:"
wc -l finches.fam
