#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p main supplement

# main figures
python scripts/fig_ibd_over_time.py
python scripts/fig_beak_shape.py
python scripts/fig_chr6_traits.py
python scripts/fig_manhattan_track.py
python scripts/fig_max_kin.py

# supplementary figures
python scripts/supp_relative_abundance.py
python scripts/supp_anc_con.py
python scripts/supp_ideogram.py
python scripts/supp_frac_founder_genome.py
python scripts/supp_kinship_top_scandens.py
python scripts/supp_ibd_length_over_time.py
python scripts/supp_drought_ibd.py
python scripts/supp_drought_ibd_bins.py
python scripts/supp_drought_beak_size.py
python scripts/supp_drought_locus_ibd.py

# drought allele frequency summary (printed, no figure)
python scripts/drought_allele_freqs.py
