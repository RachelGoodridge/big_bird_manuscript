from common import DATA, FOCAL_SPECIES, load_birds

# birds tested for each species (all six measurements and known sex), as plink --keep lists (FID 0 as in the IBIS plink files)
for spp in FOCAL_SPECIES:
    birds = load_birds(spp)
    birds.assign(FID=0)[["FID", "ind"]].to_csv(f"{DATA}/{spp}_iids.txt", sep="\t", header=False, index=False)
    print(f"{spp}: {len(birds)} birds ({(birds['sex'] == 0).sum()} male, {(birds['sex'] == 1).sum()} female)")
