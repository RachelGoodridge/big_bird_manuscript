import pandas as pd
from common import DATA, load_drought_fortis, load_master

birds = load_drought_fortis(load_master())
loci = sorted(pd.read_csv(f"{DATA}/QTLs/loci_positions.csv")["locus"].str.split("_").str[-1])

# number of large (L) alleles each bird carries at each locus, dropping predicted genotypes with < 0.99 confidence
for locus in loci:
    pred_genos = pd.read_csv(f"{DATA}/QTLs/eigenvec_df_{locus}.csv")
    pred_genos = pred_genos[pred_genos["#IID"].isin(birds["meaningful.unique"]) & (pred_genos["confidence"] >= 0.99)].set_index("#IID")
    birds[f"allele_count_{locus}"] = birds["meaningful.unique"].map(pred_genos["predicted_genotype"].str.count("L"))

# proportion of L alleles (among birds with a confident genotype) in fortis without IBD (general population) and in fortis with IBD that died
for label, group in [("without IBD", birds[~birds["has_ibd"]]), ("with IBD that died", birds[birds["has_ibd"] & (birds["outcome"] == "died")])]:
    print(f"\nfortis individuals {label} ({group.shape[0]}):")
    for locus in loci:
        if locus != "G07":  # skip G07 because it doesn't have a large allele (only blunt / pointed alleles)
            counts = group[f"allele_count_{locus}"].dropna()
            print(f"Proportion of L alleles at {locus}: {counts.sum() / (len(counts) * 2):.3f} ({len(counts)} birds genotyped)")
