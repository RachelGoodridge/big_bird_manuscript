#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p tuning tuning_plots

PLINK2=/programs/plink2_linux_20230721/plink2
IBIS=/SSD/reg259/IBIS/ibis/ibis  # v1.20.9 (commit 435a4fb)

# extract just the inds for tuning
mkdir tuning tuning_plots
$PLINK2 --bfile plink/dataset --keep dataset/birds_for_tuning.txt --autosome-num 95 --make-bed --out tuning/dataset

# create a file to store the results
echo "LENGTH,MARKERS,ERROR,POS_MEAN_ERR,NEG_MEAN_ERR,NEG_MAX_ERR" > results.csv

# loop through some of the parameter space
for LENGTH in 0.25 0.5 1 2 4; do
  for MARKERS in 50 100 200 400 800; do
    for ERROR in 0.001 0.002 0.004 0.008 0.01 0.02 0.04; do
  
      # now run IBIS
      $IBIS -b tuning/dataset -threads 16 -mL $LENGTH -mt $MARKERS -er $ERROR -printCoef -noFamID -noConvert -o tuning/finches

      # get the error rates
      ERR_RATES=$(python scripts/ibis_error.py)
      echo "$LENGTH,$MARKERS,$ERROR,$ERR_RATES" >> results.csv
      echo "finished mL = $LENGTH, mt = $MARKERS, er = $ERROR"
  
    done
  done
done

# then compare parameters
python scripts/parameter_optim.py