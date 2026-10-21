cd /SSD/reg259/big_bird_manuscript/tune_trim_ibis
FULL_VCF=/local/storage/all_birds_highcov_lowcov/autosomes_2022chrZ_all_highcov_lowcov_phased.vcf.gz
SAMPLES=Daphne_birds_for_ibis.txt
SNPS=downsampled_snps.tsv.gz
NEW_VCF=Daphne_birds_downsampled.vcf.gz
PLINK2=/programs/plink2_linux_20230721/plink2

# get all the birds from the large dataset and downsample to the SNPs in the genetic map
bcftools view --threads 16 -S $SAMPLES -T $SNPS $FULL_VCF -Oz -o $NEW_VCF
bcftools index $NEW_VCF
bcftools stats $NEW_VCF > $NEW_VCF.stats

# then use plink to filter out SNPs with >5% missing data
mkdir plink
$PLINK2 --vcf $NEW_VCF --geno 0.05 --allow-extra-chr --autosome-num 95 --make-bed --out plink/dataset

# next add the genetic map to the bim file
python make_map.py

# clear up the massive vcf file
rm $NEW_VCF $NEW_VCF.csi