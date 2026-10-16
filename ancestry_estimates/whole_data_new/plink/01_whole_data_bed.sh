#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2

# finches_pruned (built in training_new/plink/04_ld_prune.sh) was already restricted
# to exactly these 4330 samples during import, so no further --keep is needed here --
# this is a straight conversion to plink1 .bed/.bim/.fam (for ADMIXTURE) plus a PCA
# for QC/visualization. Do NOT add --geno/--maf/other frequency-dependent filters at
# this stage: the marker set must stay identical to training_new's, or ADMIXTURE's -P
# projection (which assumes matching markers in matching order) will be wrong.
$PLINK2 --pfile ../../training_new/plink/finches_pruned --allow-extra-chr --autosome-num 30 \
  --pca 10 --threads $THREADS --make-bed --out finches

echo "whole-data individuals:"
wc -l finches.fam
