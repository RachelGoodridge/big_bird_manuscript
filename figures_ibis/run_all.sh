#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p main supplement

# main figures
python scripts/fig_2_combined.py  # PC(Beak Pointedness) and IBD with the founder over time + IBD tracks along the genome
python scripts/fig_3_combined.py  # kinship to the founder's descendants + schematic of the populations (illustrations/scripts/demes_diagram.py)
python scripts/fig_4_combined.py  # beak pointedness vs IBD (genome-wide, chr6 LMM) + scandens Manhattan + chr. 6 genotype heatmap (needs ibd_phenotype_association and chr6_haplotypes output)

# supplementary figures
python scripts/supp_beak_shape.py
python scripts/supp_region_traits.py  # chr6 and chrZ regions (needs ibd_phenotype_association output)
python scripts/supp_relative_abundance.py
python scripts/supp_anc_con.py
python scripts/supp_ideogram.py
python scripts/supp_frac_founder_genome.py
python scripts/supp_founder_chrZ.py  # chrZ vs autosomal IBD with the founder for every bird with IBD
python scripts/supp_kinship_top_scandens.py
python scripts/supp_ibd_length_over_time.py
python scripts/supp_allometry.py
python scripts/supp_allometry_ground.py
python scripts/supp_morphospace.py  # body size, beak size and beak pointedness against each other
python scripts/supp_controls.py  # needs tune_trim_ibis/05_control_birds.sh
python scripts/supp_ibd_window_freqs.py
python scripts/supp_beak_traits_over_time.py
