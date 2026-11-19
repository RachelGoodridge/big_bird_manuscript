#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output

# treat the founder and each Española conirostris as a focal bird, as tune_trim_ibis/05_control_birds.sh did for the negative controls:
# segments with every other bird, keeping pairs that share more IBD than the original run's threshold (its negative control max error), then
# trimming the original run's false positive regions
python scripts/control_segments.py
