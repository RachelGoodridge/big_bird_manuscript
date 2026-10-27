import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from common import (COLOR_MAP, FOCAL_SPECIES, FOUNDER, LOCUS_IBD_COVERAGE, expand_years_alive, founder_ibd_per_bird, load_genotypes,
                    load_loci, load_master, locus_ibd_fraction, save_figure, set_plot_theme, style_legend)

set_plot_theme()

# fortis and scandens in the IBIS run, and whether each is in IBD with the founder at each locus
master = load_master()
birds = founder_ibd_per_bird(master)
birds = birds[birds["Species"].isin(FOCAL_SPECIES)].set_index("ind")
in_ibd = locus_ibd_fraction(master) >= LOCUS_IBD_COVERAGE

# one panel per locus
fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
for ax, locus in zip(axes.flat, load_loci()["name"]):
    genotypes = load_genotypes(locus)
    founder_gt = genotypes[FOUNDER]
    df = birds.assign(genotype=genotypes, in_ibd=in_ibd[locus]).dropna(subset=["genotype"])

    if founder_gt[0] == founder_gt[1]:
        # founder homozygous: every bird in IBD carries one copy of his allele from him (two only if both parents descend from him, ignored here)
        allele = founder_gt[0]
        df["allele_count"] = df["genotype"].str.count(allele)
        df["founder_copies"] = (df["in_ibd"] & (df["allele_count"] > 0)).astype(int)
        print(f"{locus} (founder {founder_gt}): {(df['in_ibd'] & (df['allele_count'] == 0)).sum()} birds in IBD without the founder's allele")
    else:
        print(f"{locus} (founder {founder_gt}): heterozygous, so only the proportion of birds in IBD is shown")

    # per year among living birds: frequency of the founder's allele (line) and of copies inherited from the founder (shaded)
    alive = expand_years_alive(df.reset_index())
    for spp in FOCAL_SPECIES:
        yearly = alive[alive["Species"] == spp].groupby("Year_alive")
        if founder_gt[0] == founder_gt[1]:
            n_alleles = yearly.size() * 2
            ax.plot(n_alleles.index, yearly["allele_count"].sum() / n_alleles, color=COLOR_MAP[spp], linewidth=1.5)
            ax.fill_between(n_alleles.index, yearly["founder_copies"].sum() / n_alleles, color=COLOR_MAP[spp], alpha=0.5, linewidth=0)
        else:
            prop = yearly["in_ibd"].mean()
            ax.fill_between(prop.index, prop.astype(float), color=COLOR_MAP[spp], alpha=0.5, linewidth=0)

    ax.set_ylabel(f"Frequency of '{allele}' Allele" if founder_gt[0] == founder_gt[1] else "Proportion of Individuals in IBD")
    ax.set_title(f"Locus {locus} (founder {founder_gt})")
    ax.set_xticks(range(1985, 2025, 10))
    ax.set_xlim(1983, 2024)  # few genotyped birds were alive before 1983
    ax.set_ylim(0, 1)
    sns.despine(top=True, right=True, ax=ax)

# one legend for all panels above them: each species' line (allele frequency) and shading (inherited from the founder)
handles = [Line2D([0], [0], color=COLOR_MAP[spp], linewidth=1.5, label=spp) for spp in FOCAL_SPECIES]
handles += [Patch(facecolor=COLOR_MAP[spp], alpha=0.5, label=rf"$\it{{{spp}}}$ inherited from the founder (IBD)") for spp in FOCAL_SPECIES]
style_legend(fig, handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False)
plt.tight_layout()
save_figure("alleles_all_loci.png", fig)
