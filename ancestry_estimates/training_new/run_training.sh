#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

# runs the full training pipeline start to finish: shared LD-pruning of the
# whole VCF, then training-specific relatedness filtering + PCA + ADMIXTURE
# across k = 2..12, plus the whole-data bed/PCA (doesn't depend on K, so no
# need to wait). Intended to be kicked off before bed and checked on in the
# morning, e.g.:
#   ./run_training.sh > run_training.log 2>&1

plink/01_subset_and_import.sh
plink/02_missingness_filter.sh
plink/03_maf_filter.sh
plink/04_ld_prune.sh
plink/05_training_bed.sh
(cd plink && Rscript 06_plot_pca.R)
admixture/01_run_admixture.sh
(cd admixture && Rscript 02_plot_admixture.R)

echo "training done. inspect admixture/finches.cv.error and (after plotting) admixture/subset_admix.png to pick k."

# also build the bed + pca for the full dataset (doesn't need K)
../whole_data_new/plink/01_whole_data_bed.sh
(cd ../whole_data_new/plink && Rscript 02_plot_pca.R)

echo "all done. once you've picked K, run whole_data_new/admixture/01_run_projection.sh <K>."
