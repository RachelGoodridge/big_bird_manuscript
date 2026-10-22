#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output tuning_plots

IBIS=/SSD/reg259/IBIS/ibis/ibis  # v1.20.9 (commit 435a4fb)

# choose -mL=0.5, -mt=200, -er=0.008
$IBIS -b plink/dataset -threads 16 -mL 0.5 -mt 200 -er 0.008 -printCoef -noFamID -noConvert -o output/finches

# have a look at the founder's known offspring
python scripts/plot_ibis.py