IBIS=/SSD/reg259/IBIS/ibis/ibis

# choose -mL=0.5, -mt=200, -er=0.008
mkdir output
$IBIS -b plink/dataset -threads 16 -mL 0.5 -mt 200 -er 0.008 -printCoef -noFamID -noConvert -o output/finches

# have a look at the founder's known offspring
python plot_ibis.py