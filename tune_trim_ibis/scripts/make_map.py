import os
from os.path import join
import pandas as pd
from pandas_plink import read_plink
import warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="pandera")

# get the plink files
(bim, fam, bed) = read_plink("plink/dataset*.bed", verbose=False)

# update the chromosome names
chromosomes = list(bim["chrom"].unique())
chr_map = {chrom: int(i+1) for i, chrom in enumerate(chromosomes)}
bim["chrom"] = bim["chrom"].map(chr_map)

# add snp identifiers
bim["snp"] = bim["chrom"].astype(str) + ":" + bim["pos"].astype(str)

# read in the gmap files
gmap_path = "/local/storage/gmap_files"
gmap_list = []
for gmap_file in os.listdir(gmap_path):
    
    # figure out the chromosome number
    if gmap_file[3:-5] in chr_map:
        chrom = chr_map[gmap_file[3:-5]]
    elif gmap_file[:-5] in chr_map:
        chrom = chr_map[gmap_file[:-5]]
    else:
        print(f"Chromosome not found for gmap file: {gmap_file}")
        continue
    
    # add the genetic map to a list
    data = pd.read_csv(join(gmap_path, gmap_file), sep="\t")
    data["snp"] = str(chrom) + ":" + data["pos"].astype(str)
    gmap_list.append(data)

# add the genetic map to the bim file
gmap_data = pd.concat(gmap_list)
bim["cm"] = bim["snp"].map(gmap_data.set_index("snp")["cM"])

# output the new bim file
assert not bim.isna().any().any()
bim.drop(bim.columns[-1], axis=1, inplace=True)
bim.to_csv("plink/dataset.bim", sep="\t", index=False, header=False)

# also save the chromosome map as a text file
with open("dataset/chr_map.txt", "w") as f:
    for key, value in chr_map.items():
        f.write(f"{key} {value}\n")