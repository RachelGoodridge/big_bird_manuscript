#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# hand-made panels, each also saved on its own
python make_galapagos_map.py   # galapagos_map.svg (figure 1A; needs galapagos_land.geojson and conirostris_head.png)
python pedigree.py             # figure1_pedigree.svg (figure 1B)
python timeline.py             # figure1_timeline.svg (alternative figure 1B)
python demes_diagram.py        # figure3_demes_diagram.svg (schematic for figure 3)

# assembled figures
python fig_1_combined.py       # figure1.svg (map + pedigree + ADMIXTURE)
python fig_1v2_combined.py     # figure1v2.svg (map + timeline + ADMIXTURE)
