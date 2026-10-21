#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

# set this to whatever K you pick from training_new/admixture/finches.cv.error +
# subset_admix.png (and cv_error.png)
K=9

THREADS=32
ADMIXTURE=/programs/admixture/admixture
TRAINING_ADMIX=../../training_new/admixture

cp ../plink/finches.bed ../plink/finches.bim ../plink/finches.fam .

# ADMIXTURE doesn't accept non-human chromosome codes -- zero out the first column
awk '{$1="0";print $0}' finches.bim > finches.bim.tmp
mv finches.bim.tmp finches.bim

# lock the study set's allele frequencies to the chosen training run's P matrix and
# project every individual onto it (this is why the marker sets had to match exactly)
cp $TRAINING_ADMIX/finches.${K}.P finches.${K}.P.in
$ADMIXTURE -P finches.bed $K -j$THREADS > log${K}.out

echo "done -- ancestry proportions in finches.${K}.Q"
