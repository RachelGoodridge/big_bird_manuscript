#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data output figures

# needs the trimmed IBD segments and plink files from tune_trim_ibis
./01_make_grm.sh  # birds tested and their LOCO GRMs (slow; only needs rerunning if the birds or SNPs change)
python scripts/run_association.py
python scripts/plot_manhattan.py
python scripts/plot_diagnostics.py
