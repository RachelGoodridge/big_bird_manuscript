#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p figures

# needs the trimmed IBD segments from tune_trim_ibis and the genotype calls from locus_genotyping
python scripts/drought_boxplots.py
python scripts/drought_ibd_bins.py
python scripts/drought_locus_ibd.py

# drought allele frequency summary (printed, no figure)
python scripts/drought_allele_freqs.py
