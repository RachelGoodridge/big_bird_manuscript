#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p dataset

FULL_VCF=/local/storage/all_birds_highcov_lowcov/autosomes_2022chrZ_all_highcov_lowcov_phased.vcf.gz
GMAP_DIR=/local/storage/gmap_files
SNPS=dataset/downsampled_snps.tsv.gz

# work in a temporary directory that is removed on exit
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# get every position in the genetic maps (columns: pos, chr, cM)
for GMAP in $GMAP_DIR/*.gmap; do
  tail -n +2 $GMAP | awk -v OFS='\t' '{print $2, $1}'
done | sort -k1,1 -k2,2n > $TMP/gmap_pos.tsv

# keep only the genetic map SNPs that are present in the full vcf
bcftools query -T $TMP/gmap_pos.tsv -f '%CHROM\t%POS\n' $FULL_VCF | bgzip > $TMP/snps.tsv.gz

# never replace an existing SNP list with a different one
if [[ -e $SNPS ]]; then
  if cmp -s <(zcat $SNPS | sort) <(zcat $TMP/snps.tsv.gz | sort); then
    echo "$SNPS already lists the same $(zcat $SNPS | wc -l) SNPs, keeping it unchanged"
    exit 0
  fi
  mv $TMP/snps.tsv.gz $SNPS.new
  echo "ERROR: the new SNP list differs from $SNPS, wrote it to $SNPS.new instead" >&2
  exit 1
fi
mv $TMP/snps.tsv.gz $SNPS
tabix -s1 -b2 -e2 $SNPS
