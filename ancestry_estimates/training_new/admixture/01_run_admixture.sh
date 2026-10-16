#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

THREADS=32
ADMIXTURE=/programs/admixture/admixture

cp ../plink/finches.bed ../plink/finches.bim ../plink/finches.fam .

# ADMIXTURE doesn't accept non-human chromosome codes -- zero out the first column
awk '{$1="0";print $0}' finches.bim > finches.bim.tmp
mv finches.bim.tmp finches.bim

for k in {2..12}; do
  echo "running k=$k..."
  $ADMIXTURE --cv finches.bed $k -j$THREADS > log${k}.out
done

# collect the cross-validation error curve (lower is better, but see the note in
# plot_admixture.R about not just picking the CV-minimum for this kind of panel)
awk '/CV/ {print $3, $4}' log*.out | sed 's/(K=\([0-9]\+\)):/\1/' | sort -n -k1,1 > finches.cv.error
cat finches.cv.error
