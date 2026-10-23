#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output

# repeat the founder's filtering and trimming for each negative control (as if it were the founder): segments with every
# other bird, keeping pairs that share more IBD than the negative control max error, then trimming the false positive regions
python scripts/control_segments.py
