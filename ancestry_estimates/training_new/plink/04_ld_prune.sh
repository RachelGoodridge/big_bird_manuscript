#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2

# use standard pruning settings
$PLINK2 --pfile finches_02_maf --allow-extra-chr --autosome-num 30 \
  --indep-pairwise 50 10 0.1 --threads $THREADS --out finches_prune

# extract the pruned marker set
$PLINK2 --pfile finches_02_maf --allow-extra-chr --autosome-num 30 \
  --extract finches_prune.prune.in --threads $THREADS --make-pgen --out finches_pruned

echo "variants remaining after pruning:"
wc -l finches_pruned.pvar
