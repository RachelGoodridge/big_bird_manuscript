#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

PLINK2=/programs/plink2_linux_20230721/plink2

# the two groups: conirostris, the founder, and magnirostris vs scandens without IBD
python scripts/fst_groups.py

# per-SNP Hudson FST between the groups (all SNPs in the region, not only the pruned ones), and each group's allele counts for windowed FST
tail -n +2 data/fst_groups.tsv | cut -f1,2 > work/fst_keep.txt
$PLINK2 --bfile work/region --keep work/fst_keep.txt --pheno data/fst_groups.tsv --fst group report-variants --out work/fst --silent
for GROUP in con_mag scandens; do
  awk -v g=$GROUP '$3 == g' data/fst_groups.tsv > work/keep_$GROUP.txt
  $PLINK2 --bfile work/region --keep work/keep_$GROUP.txt --freq counts --out work/freq_$GROUP --silent
done

# FST in sliding windows along the region, and its relationship with recombination
python scripts/fst_windows.py
python scripts/fst_vs_recombination.py
