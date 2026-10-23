#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p trimming trimming_plots output

# exclude regions where more than 5% of negative controls share IBD with the founder (likely false positives),
# and trim them out of the founder's segments
python scripts/trim_segments.py

# plot the false positive regions and the locations of the trimmed segments
python scripts/plot_false_pos_loc.py
python scripts/plot_IBD_loc.py
