import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common import (COLOR_MAP, FOCAL_SPECIES, add_bird_info, expand_years_alive, format_year_axis, founder_ibd_per_bird,
                    load_bim, load_master, load_trimmed_segments, save_figure, set_plot_theme, style_legend, yearly_summary)

set_plot_theme()


# fraction of all SNPs that fall inside at least one segment (segment positions in Mb)
def fraction_snps_covered(segs, bim):
    covered = 0
    for chrom, pos in bim.groupby("chrom")["pos"]:
        pos = np.sort(pos.to_numpy())
        chrom_segs = segs[segs["chrom"] == chrom]
        depth = np.zeros(len(pos) + 1, dtype=int)
        np.add.at(depth, np.searchsorted(pos, np.rint(chrom_segs["phys_start_pos"] * 1e6), side="left"), 1)
        np.add.at(depth, np.searchsorted(pos, np.rint(chrom_segs["phys_end_pos"] * 1e6), side="right"), -1)
        covered += (np.cumsum(depth)[:-1] > 0).sum()
    return covered / len(bim)


########################### LOAD DATA ###########################
master = load_master()
bim = load_bim()

# trimmed segments of each bird for every year it was alive
seg_alive = expand_years_alive(add_bird_info(load_trimmed_segments(), master))
years = sorted(seg_alive["Year_alive"].unique())
seg_alive = seg_alive[seg_alive["Species"].isin(FOCAL_SPECIES)]

# fraction of the founder's genome (SNPs) carried by at least one living fortis / scandens each year
frac_genome = pd.DataFrame([{"Year": year, "Species": spp, "frac": fraction_snps_covered(seg_alive[(seg_alive["Year_alive"] == year) & (seg_alive["Species"] == spp)], bim)}
                            for year in years for spp in FOCAL_SPECIES])

# birds alive (count_all) and birds carrying IBD (count_ibd) each year, for the point sizes
frac_genome = pd.merge(frac_genome, yearly_summary(founder_ibd_per_bird(master)), on=["Year", "Species"], how="left").fillna(0)


########################### CREATE PLOT ###########################
fig, ax = plt.subplots(1, 1, figsize=(6, 4.875))
for spp in FOCAL_SPECIES:
    sub = frac_genome[frac_genome["Species"] == spp]
    sns.lineplot(x="Year", y="frac", data=sub, label=spp, color=COLOR_MAP[spp], ax=ax, linewidth=1.5, zorder=2)

    # filled points scale with the number of carriers, open circles with the number of birds alive
    ax.scatter(sub["Year"], sub["frac"], s=sub["count_ibd"] * 5, color=COLOR_MAP[spp], alpha=0.8, linewidth=0)
    ax.scatter(sub["Year"], sub["frac"], s=sub["count_all"] * 5, facecolors="none", edgecolors=COLOR_MAP[spp], linewidth=0.8)

ax.set_ylabel("Fraction of the Founder's Genome Retained")
format_year_axis(ax)
ax.set_xlim(1983, 2024)
ax.set_ylim(-0.12, 1.05)
style_legend(ax, loc="lower right")
sns.despine(top=True, right=True, left=False, bottom=False, ax=ax)

plt.tight_layout()
save_figure("supplement/frac_founder_exist.png")
