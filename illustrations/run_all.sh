#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/scripts"

# hand-made panels, each also saved on its own in figures/ (inputs in data/)
python make_galapagos_map.py   # galapagos_map.svg (figure 1A; needs data/galapagos_land.geojson and data/conirostris_head.png)
python pedigree.py             # founder_pedigree.svg (simplified pedigree of the Big Bird lineage, on its own)
python timeline.py             # founder_timeline.svg (figure 1B)
python demes_diagram.py        # demes_diagram.svg (schematic for figure 3)

# assembled figure
python fig_1_combined.py       # figure1.svg / figure1.pdf (map + timeline + ADMIXTURE)
