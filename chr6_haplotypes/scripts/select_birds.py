import pandas as pd
from common import (CANDIDATES, CHROM, DATA, END, GROUPS, MIN_ANCESTRY, REGION_COVERAGE, START, VCF_SAMPLES, load_founder_coef, load_master, region_ibd,
                    species_ancestry)

# candidate birds of each group (before the regional check, regional_check.py)

# birds in the phased VCF, with the newest lineage classifications
master = load_master(island=None)
in_vcf = set(pd.read_csv(VCF_SAMPLES, header=None)[0])
master = master[master["meaningful.unique"].isin(in_vcf)].rename(columns={"meaningful.unique": "ind"})

# IBD with the founder in the chromosome 6 region for the Daphne fortis and scandens (all are in the IBIS run)
daphne = master[(master["Island"] == "Daphne") & master["Species"].isin(["fortis", "scandens"])].copy()
assert daphne["ind"].isin(load_founder_coef()["ind"]).all(), "some Daphne fortis / scandens are not in the IBIS run"
daphne["region_ibd"] = region_ibd(daphne["ind"].to_numpy())
ibd_rule = {
    "fortis": daphne["Species"].eq("fortis") & daphne["region_ibd"].eq(0),                                # fortis without any IBD in the region
    "scandens_no_ibd": daphne["Species"].eq("scandens") & daphne["region_ibd"].eq(0),                     # scandens without any
    "scandens_ibd": daphne["Species"].eq("scandens") & daphne["region_ibd"].ge(REGION_COVERAGE),          # scandens with >= 50% (as in figure 4B)
}

# each group's birds with at least MIN_ANCESTRY genome-wide ancestry of their own species (magnirostris: both magnirostris clusters summed), except the
# scandens with IBD, who carry Big Bird ancestry across the genome by definition (all kept, as in figure 4B); then those with the most ancestry, up
# to the group's maximum
ancestry = {spp: species_ancestry(spp) for spp in GROUPS["ancestry"].unique()}
birds = []
for group, g in GROUPS.iterrows():
    if group in ibd_rule:
        candidates = daphne[ibd_rule[group]]
    else:
        candidates = master[(master["Island"] == g["island"]) & (master["Species"] == g["species"])]
    candidates = candidates.sort_values("ind").assign(ancestry=lambda d: d["ind"].map(ancestry[g["ancestry"]]))
    n = len(candidates)
    low = candidates[(candidates["ancestry"] < MIN_ANCESTRY) & (group != "scandens_ibd")]
    candidates = candidates.drop(low.index)
    if pd.notna(g["max_birds"]):
        candidates = candidates.nlargest(int(g["max_birds"]), "ancestry").sort_values("ind")
    print(f"{group}: {len(candidates)} of {n} birds ({len(low)} with {g['ancestry']} ancestry < {MIN_ANCESTRY} left out), ancestry of those kept "
          f"{candidates['ancestry'].min():.2f}-{candidates['ancestry'].max():.2f}")
    birds.append(pd.DataFrame({"ind": candidates["ind"], "group": group, "ancestry": candidates["ancestry"]}))
birds = pd.concat(birds, ignore_index=True)
birds.insert(2, "region_ibd", birds["ind"].map(daphne.set_index("ind")["region_ibd"]))  # Daphne fortis / scandens only

# partial carriers left out (as in the region test)
partial = daphne[daphne["region_ibd"].between(0, REGION_COVERAGE, inclusive="neither")]
print(f"left out: {(partial['Species'] == 'scandens').sum()} scandens with 0-{REGION_COVERAGE:.0%} of chr{CHROM}:{START}-{END} in IBD, "
      f"{(daphne['Species'].eq('fortis') & daphne['region_ibd'].gt(0)).sum()} fortis with any")

birds.to_csv(CANDIDATES, sep="\t", index=False)
birds["ind"].to_csv(f"{DATA}/candidates.txt", index=False, header=False)
print(f"saved {len(birds)} candidate birds to {CANDIDATES}")
