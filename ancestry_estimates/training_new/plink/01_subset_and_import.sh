#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

# ── shared preprocessing (steps 01-04) ──────────────────────────────────────
# Everything through 04_ld_prune.sh operates on the union of every sample either
# project will ever use: whole_data_new/plink/bird_selection.txt (4330 samples,
# which is training's keep_list.txt (156) union every Daphne bird in the master
# database). Restricting to this set *before* computing missingness/MAF/LD
# statistics keeps those statistics from being skewed by the ~468 other samples
# in the raw VCF that neither project touches (other unrelated taxa entirely).
# The resulting pruned marker set (finches_pruned.pgen) is the canonical SNP
# panel both training_new/plink/05_training_bed.sh and whole_data_new build
# their .bed files from -- do not re-filter by allele frequency downstream of
# that point, or the training and whole-data marker sets will diverge, which
# breaks ADMIXTURE's -P projection (it requires identical markers in identical
# order between the reference .P file and the study .bed/.bim).

THREADS=32
PLINK2=/programs/plink2_linux_20230721/plink2
FULL_VCF=/local/storage/all_birds_highcov_lowcov/autosomes_2022chrZ_all_highcov_lowcov_phased.vcf.gz

# chrZ is excluded: it's hemizygous in females (ZW system), and if female
# genotype calls were forced diploid at those sites that introduces artificial
# homozygosity which would bias kinship/PCA/admixture. Matches the old
# autosomes-only workflow.
$PLINK2 --vcf $FULL_VCF --not-chr chrZ --allow-extra-chr --autosome-num 30 \
  --set-missing-var-ids @:# --threads $THREADS --make-pgen --out finches_00_raw
