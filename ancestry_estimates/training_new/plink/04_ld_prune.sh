#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2

# identify a quasi-independent marker set: 50kb window, step 10 variants, r^2 < 0.1.
# (the exact parameters used in the previous round of this analysis were lost, so
# this is a standard default for ADMIXTURE/PCA-style pruning -- adjust if desired)
$PLINK2 --pfile finches_02_maf --allow-extra-chr --autosome-num 30 \
  --indep-pairwise 50 10 0.1 --threads $THREADS --out finches_prune

# extract the pruned marker set. this is the canonical SNP panel + sample set
# (4330 individuals) that both training_new/plink/05_training_bed.sh and
# whole_data_new/plink/01_whole_data_bed.sh build their .bed files from.
$PLINK2 --pfile finches_02_maf --allow-extra-chr --autosome-num 30 \
  --extract finches_prune.prune.in --threads $THREADS --make-pgen --out finches_pruned

echo "variants remaining after pruning:"
wc -l finches_pruned.pvar
