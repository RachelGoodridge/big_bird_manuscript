import os
import sys
import pandas as pd
import shutil

# move to the directory of this script
os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))

# get a list of the reference panel
panel="/local/storage/unphased_reference_panel/finches_sentieon_concat_SNPs_fil_setGT_PASS_BIALLELIC_ACf.csv"
data = pd.read_csv(panel)

# check out the birds from big bird paper admixture plot (Fig 2C)
spp_island = {"difficilis":["Fernandina","Santiago","Pinta"], "septentrionalis":["Wenman(Wolf)","Culpepper(Darwin)","Darwin"], "conirostris":["Espanola"], "propinqua":["Genovesa"],
              "acutirostris":["Genovesa"], "fortis":["Daphne","Santa_Cruz"], "fuliginosa":["Santiago","Santa_Cruz"], "magnirostris":["Genovesa", "Daphne"], "scandens":["Daphne"]}

# write selected birds to a file
file_path="keep_list.txt"
for key, value in spp_island.items():
    mask = (data["Species_inferred"]==key) & (data["Island_inferred"].isin(value))
    with open(file_path, "a") as file:
        for bird in data[mask]["meaningful.unique"].to_list():
            file.write(bird + "\n")

# get a list of the older Daphne birds
old="/local/storage/Daphne_pruned/Daphne_neutral_aims_pruned_not2022.csv"
data = pd.read_csv(old)

# select the 10 of each species with highest self-ancestry
species = ["fortis", "magnirostris", "scandens", "fuliginosa", "magnirostris"]
ancestries = ["anc_fort", "anc_mag", "anc_scan", "anc_fuli", "anc_mag2"]
for spp, anc in zip(species, ancestries):
    birds = data[data["Species"] == spp].sort_values(by=anc, ascending=False).head(10)["sample"].to_list()
    with open(file_path, "a") as file:
        for bird in birds:
            file.write(bird + "\n")

# add the founder to the list
with open(file_path, "a") as file:
    file.write("04Dap5110\n")

# remove the duplicates
with open(file_path, "r") as file:
    lines = file.readlines()
lines = sorted(list(set(lines)))
with open(file_path, "w") as file:
    for line in lines:
        file.write(line)