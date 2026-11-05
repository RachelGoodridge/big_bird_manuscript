#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data grm

PLINK2=/programs/plink2_linux_20230721/plink2
GEMMA=/programs/gemma_0.98.5/gemma  # v0.98.5
DATASET=../tune_trim_ibis/plink/dataset  # the IBIS SNPs (genotyped in >= 95% of birds)
N_PARALLEL=15  # GRMs built at once
export OPENBLAS_NUM_THREADS=4

# birds tested for each species (all six measurements and known sex)
python scripts/select_birds.py

for SPP in fortis scandens; do
  rm -rf grm/$SPP
  mkdir -p grm/$SPP

  # genotypes of just this species
  $PLINK2 --bfile $DATASET --keep data/${SPP}_iids.txt --autosome-num 95 --make-bed --out grm/$SPP/plink

  # GEMMA needs a phenotype to build a relatedness matrix, so give every bird the same dummy value
  awk '{print 1}' grm/$SPP/plink.fam > grm/$SPP/dummy_pheno.txt

  # one leave-one-chromosome-out (LOCO) centered relatedness matrix per chromosome, from the SNPs on all other chromosomes
  CHROMS=$(cut -f1 grm/$SPP/plink.bim | sort -un)
  for CHR in $CHROMS; do
    awk -v chr="$CHR" '$1 != chr {print $2}' grm/$SPP/plink.bim > grm/$SPP/snps_excl_chr$CHR.txt
  done
  echo $CHROMS | tr ' ' '\n' | xargs -P $N_PARALLEL -I{} sh -c \
    "$GEMMA -bfile grm/$SPP/plink -p grm/$SPP/dummy_pheno.txt -gk 1 -snps grm/$SPP/snps_excl_chr{}.txt -outdir grm/$SPP -o loco_chr{} > /dev/null && echo '$SPP chr{} done'"
  rm grm/$SPP/snps_excl_chr*.txt
done
