import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from common import (COLOR_MAP, GROUPS, MAX_CHRZ_MISSING, chrz_missing_birds, founder_ibd_per_bird, ibd_label, load_bim, load_chromosomes,
                    load_master, load_trimmed_segments, order_groups, save_figure, set_plot_theme, species_label, style_legend)

# does the founder's Z chromosome persist in his descendants? Each bird's fraction of chrZ SNPs in (trimmed) IBD with the founder against its
# fraction of autosomal SNPs, for every bird with any IBD. Birds that inherited the founder's Z along with the rest of his genome fall near the
# diagonal; a Z inherited without much else (or lost while the autosomes were kept) moves a bird above (or below) it. Males (two Z) and females
# (one Z, always from their father) are told apart, as a female can only get the founder's Z through her father
set_plot_theme()
master = load_master()
chrz = next(code for code, name in load_chromosomes().items() if name == "Z")

# SNPs on the autosomes and on chrZ, and each bird's SNPs in IBD on each
bim = load_bim()
n_snps = {"auto": (bim["chrom"] != chrz).sum(), "z": (bim["chrom"] == chrz).sum()}
seg = load_trimmed_segments()
ibd = {"auto": seg[seg["chrom"] != chrz].groupby("ind")["SNP_count"].sum(), "z": seg[seg["chrom"] == chrz].groupby("ind")["SNP_count"].sum()}

birds = founder_ibd_per_bird(master)
birds = order_groups(birds[birds["fraction_snps_ibd"] > 0], GROUPS)
for part in ["auto", "z"]:
    birds[f"ibd_{part}"] = (birds["ind"].map(ibd[part]) / n_snps[part]).fillna(0)
birds["sex"] = birds["ind"].map(master.set_index("meaningful.unique")["sex"]).map({0: "male", 1: "female"}).fillna("unknown")

# birds with much of chrZ missing are left out, as in the association test (MAX_CHRZ_MISSING, common.py)
missing = birds["ind"].isin(chrz_missing_birds())
print(f"left out (more than {MAX_CHRZ_MISSING:.0%} of chrZ missing): {missing.sum()} birds, {birds.loc[missing, 'Species'].value_counts().to_dict()}")
birds = birds[~missing]

# summary: birds of each group with the founder's Z, and how much of it
summary = birds.groupby(["Species", "sex"], observed=True).agg(birds=("ind", "size"), with_z=("ibd_z", lambda z: (z > 0).sum()),
                                                                mean_z=("ibd_z", "mean"), mean_auto=("ibd_auto", "mean"))
print(summary.round(3).to_string())
recent = birds[(birds["Species"] == "scandens") & (birds["Year"] >= 2011)].nlargest(20, "ibd_auto")
print(f"\nthe 20 scandens first seen 2011-2022 with the most autosomal IBD: {(recent['ibd_z'] > 0).sum()} carry some of the founder's Z "
      f"(mean fraction {recent['ibd_z'].mean():.2f}; autosomal {recent['ibd_auto'].mean():.2f})")


########################### CREATE PLOT ###########################
MARKERS = {"female": "o", "male": "s", "unknown": "^"}  # as in a pedigree, with triangles for birds of unknown sex
fig, ax = plt.subplots(1, 1, figsize=(9.3, 6.5))
ax.set_box_aspect(1)  # square plot area, with the legend to its right
# fortis and scandens first, so the founder's lineages are drawn on top of them
for group in GROUPS:
    for sex, marker in MARKERS.items():
        sub = birds[(birds["Species"] == group) & (birds["sex"] == sex)]
        ax.scatter(sub["ibd_auto"], sub["ibd_z"], color=COLOR_MAP[group], marker=marker, s=28, alpha=0.7, linewidth=0)
ax.plot([0, 1], [0, 1], color="gray", linestyle="--", linewidth=1, zorder=0)
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.set_xlabel(ibd_label("Proportion of autosomal SNPs in IBD"))
ax.set_ylabel(ibd_label("Proportion of chrZ SNPs in IBD"))

handles = [Line2D([], [], marker="o", color="w", markerfacecolor=COLOR_MAP[group], markersize=7, label=species_label(group)) for group in GROUPS]
handles += [Line2D([], [], marker=MARKERS[sex], color="w", markerfacecolor="gray", markersize=7, label=sex) for sex in MARKERS]
style_legend(ax, handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5))
sns.despine(top=True, right=True, left=False, bottom=False)
plt.tight_layout()
save_figure("supplement/founder_chrZ_vs_autosomes.png")
