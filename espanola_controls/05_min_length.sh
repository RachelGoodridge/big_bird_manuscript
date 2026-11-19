#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# test of a minimum IBD segment length, applied after the IBD threshold and before trimming the false positive regions: each pair's segments
# separated by gaps of at most 1 cM are merged (IBIS often splits one long segment into short pieces), and segments must be in merged blocks
# longer than the longest one the founder shares with any negative control (5.375 cM). Outputs in min_length_merge1cM/
mkdir -p min_length_merge1cM
python scripts/control_segments.py --min-length
python scripts/plots.py --min-length
python scripts/plot_tracks.py --min-length
