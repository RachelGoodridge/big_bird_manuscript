#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output

IBIS=/SSD/reg259/ibis/ibis  # v1.20.9 (commit 435a4fb), as in the original run

# the original run's parameters (tune_trim_ibis/03_run_ibis.sh: -mL 0.5, -mt 200, -er 0.008), comparing only the first samples (the founder
# and the Española conirostris, indices 0 to N-1) with every bird
LAST=$(wc -l < dataset/focal_birds.txt)  # -setIndexEnd is exclusive
$IBIS -b plink/dataset -threads 16 -mL 0.5 -mt 200 -er 0.008 -printCoef -noFamID -noConvert -setIndexStart 0 -setIndexEnd $LAST -o output/espanola

# check that the founder's segments and kinship are exactly those of the original run (so the two runs are comparable)
python scripts/check_founder.py
