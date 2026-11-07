import pandas as pd
from common import DONORS, FST_GROUPS, load_birds

# the two FST groups (plink --pheno file): the donors (conirostris, the founder, and magnirostris) vs the scandens without IBD
birds = load_birds()
fst_groups = pd.concat([pd.DataFrame({"#FID": 0, "IID": birds.loc[birds["group"].isin(DONORS), "ind"], "group": "con_mag"}),
                        pd.DataFrame({"#FID": 0, "IID": birds.loc[birds["group"] == "scandens_no_ibd", "ind"], "group": "scandens"})])
fst_groups.to_csv(FST_GROUPS, sep="\t", index=False)
print(fst_groups["group"].value_counts().to_dict())
