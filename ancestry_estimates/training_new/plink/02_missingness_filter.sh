#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2

# drop variants missing in > 10% of samples; adjust --geno if you want stricter/looser.
# (no --mind sample-missingness filter here on purpose -- both bird_selection.txt and
# keep_list.txt are hand-curated sample lists, and we don't want this step to silently
# drop one of those samples for being low-coverage)
$PLINK2 --pfile finches_00_raw --allow-extra-chr --autosome-num 30 \
  --geno 0.1 --threads $THREADS --make-pgen --out finches_01_geno
