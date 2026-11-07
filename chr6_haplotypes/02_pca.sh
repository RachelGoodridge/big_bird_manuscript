#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

PLINK2=/programs/plink2_linux_20230721/plink2

# LD-pruned SNPs (r2 < 0.2 in 50 SNP windows, step 10) of the birds analyzed and their PCA (genotypes)
$PLINK2 --bfile work/region --indep-pairwise 50 10 0.2 --out work/pruned --silent
$PLINK2 --bfile work/region --extract work/pruned.prune.in --pca 10 --out work/pca --silent
echo "$(wc -l < work/pruned.prune.in) of $(wc -l < work/region.bim) SNPs left after LD pruning"

# PCA colored by group, and each group's mean position
python scripts/plot_pca.py
