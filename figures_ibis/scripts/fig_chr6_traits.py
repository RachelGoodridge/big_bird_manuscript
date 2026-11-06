import importlib.util
import sys
import pandas as pd
import matplotlib.pyplot as plt
from common import ASSOCIATION, load_chromosomes, save_figure, set_plot_theme

# the IBD-phenotype association helpers (loaded by path since it is also named common.py)
_spec = importlib.util.spec_from_file_location("association_common", f"{ASSOCIATION}/scripts/common.py")
sys.modules["association_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["association_common"])
from association_common import CHR6_REGION, TRAITS, plot_region_box, region_lmm, region_test_birds

set_plot_theme()
spp = "scandens"
chrom_name, start, end = CHR6_REGION

# the scandens of the IBD-phenotype association test (traits fit within them) with at least half of the chromosome 6 region in IBD, or none
chrom_code = {name: code for code, name in load_chromosomes().items()}[chrom_name]
birds, eig = region_test_birds(spp, chrom_code, start, end)

# LMM of each trait on the IBD group, with the same LOCO GRM and covariates (sex, IBD outside the chromosome, and PC(Body Size) for beak size /
# pointedness) as the windows of the IBD-phenotype association
stats = pd.DataFrame({trait: region_lmm(birds, eig, trait) for trait in TRAITS}, index=["beta", "p_value"]).T
print(birds["ibd_group"].value_counts().to_dict())
print(stats)

fig, axes = plt.subplots(2, 2, figsize=(10.6, 10))
for ax, trait in zip(axes.flat, TRAITS):
    plot_region_box(ax, birds, trait, stats.loc[trait, "beta"], stats.loc[trait, "p_value"], spp)
fig.suptitle(f"IBD on chr. {chrom_name} ({start / 1e6:g}-{end / 1e6:g} Mb) in " + r"$\it{" + spp + "}$")
plt.tight_layout(h_pad=3, w_pad=3)  # extra room between panels
save_figure(f"main/chr{chrom_name}_IBD_vs_traits.png")
