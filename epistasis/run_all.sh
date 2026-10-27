#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p figures

# needs the genotype calls from locus_genotyping/run_all.sh
python scripts/plot_genotype_by_ibd.py
python scripts/plot_fortis_interaction.py
