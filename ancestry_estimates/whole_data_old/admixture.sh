# change to correct directory
cd /SSD/reg259/expanded_big_bird_analysis/Daphne_project

# define files needed for the training subset
REF_10=/SSD/reg259/expanded_big_bird_analysis/training_subset/admixture/finches.10.P
REF_11=/SSD/reg259/expanded_big_bird_analysis/training_subset/admixture/finches.11.P

# define files needed for the whole data
STUDY_BIM=/SSD/reg259/expanded_big_bird_analysis/Daphne_plink/finches.bim
STUDY_BED=/SSD/reg259/expanded_big_bird_analysis/Daphne_plink/finches.bed
STUDY_FAM=/SSD/reg259/expanded_big_bird_analysis/Daphne_plink/finches.fam

# copy files to the new folder
cp $STUDY_BIM finches.bim
cp $STUDY_BED finches.bed
cp $STUDY_FAM finches.fam

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' finches.bim > finches.bim.tmp
mv finches.bim.tmp finches.bim

# run projection for k=10
cp $REF_10 finches.10.P.in
/programs/admixture/admixture -P finches.bed 10 -j20 > log10.out

# run projection for k=11
cp $REF_11 finches.11.P.in
/programs/admixture/admixture -P finches.bed 11 -j20 > log11.out