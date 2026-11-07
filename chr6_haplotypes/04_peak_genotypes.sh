#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# genotypes at the diagnostic SNPs in the FST peaks (output/fst_peak_regions.tsv), as copies of the conirostris / magnirostris allele
python scripts/plot_peak_heatmap.py
