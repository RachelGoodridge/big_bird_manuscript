#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p plots

# (1) like figure 2B: each Daphne bird's IBD with the founder and with each Española conirostris
python scripts/plots.py
# (2) like figure 2E: where along the genome each focal bird shares IBD with Daphne fortis / scandens
python scripts/plot_tracks.py
