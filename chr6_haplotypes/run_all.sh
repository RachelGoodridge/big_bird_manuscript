#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# needs the trimmed IBD segments (tune_trim_ibis) and the chromosome 6 region (CHR6_REGION, ibd_phenotype_association/scripts/common.py)
./01_extract_region.sh   # birds compared (IBD, ancestry and regional checks; fortis only as a reference), and the chromosome 6 region for them (reads the full VCF)
./02_pca.sh              # PCA of the LD-pruned SNPs
./03_fst.sh              # FST between conirostris + magnirostris and scandens without IBD, and its relationship with recombination
./04_peak_genotypes.sh   # genotypes at the diagnostic SNPs in the FST peaks
