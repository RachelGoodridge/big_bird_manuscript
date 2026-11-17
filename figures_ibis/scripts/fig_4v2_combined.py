# alternative to figure 4 (fig_4_combined.py) with the chromosome 6 genotype heatmap (chr6_haplotypes/figures/peak_heatmap.png) as panel D:
# A, C (narrower) and the heatmap (D) in the top row, the Manhattan plot (B) below them
import importlib.util
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
from common import (ASSOCIATION, add_panel_letters, load_chromosomes, load_ibd_beak_birds, load_master, plot_beak_vs_ibd, save_figure, set_plot_theme,
                    species_label, style_legend)

# the IBD-phenotype association helpers (loaded by path since it is also named common.py)
_spec = importlib.util.spec_from_file_location("association_common", f"{ASSOCIATION}/scripts/common.py")
sys.modules["association_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["association_common"])
from association_common import CHR6_REGION, genome_layout, load_results, plot_manhattan, plot_region_box, region_lmm, region_test_birds

from ibis_common import IBD_F

# the chromosome 6 genotype heatmap (chr6_haplotypes), loaded by path for the same reason
_spec = importlib.util.spec_from_file_location("haplotype_common", f"{os.path.dirname(ASSOCIATION)}/chr6_haplotypes/scripts/common.py")
sys.modules["haplotype_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["haplotype_common"])
from haplotype_common import GROUPS, HEATMAP_SNPS, PEAK_REGIONS, plot_peak_heatmap

set_plot_theme()
SPP = "scandens"
TRAIT = "PC(Beak Pointedness)"
CHROM, START, END = CHR6_REGION


# layout, in inches. Panels A and D keep the proportions they have in figure 4 and in the standalone heatmap
ROW_HEIGHT = 4.6             # top row
A_WIDTH = ROW_HEIGHT * 1.21  # as in figure 4
C_WIDTH = 2.5
TREE_WIDTH, GROUPS_WIDTH, HEATMAP_WIDTH = 0.58, 0.12, ROW_HEIGHT * 1.04  # the heatmap's tree, group colors and genotypes, as in the standalone plot
LEGEND_WIDTH = 1.95          # room right of the heatmap for its legends
GAP_AC, GAP_CD = 1.3, 0.9  # between the panels of the top row
B_ASPECT = 5.0               # the Manhattan plot's width / height (3.63 in figure 4)
LEFT, RIGHT, TOP, BOTTOM, ROW_GAP = 1.0, 0.1, 1.0, 0.8, 1.3  # margins, and the gap between the rows
HEATMAP_FONTSIZE = 13
HEATMAP_COLORS = ["#ececec", "#9E9E9E", "#505050"]  # 0, 1, 2 copies of the donor allele: grays, as the groups beside the heatmap are in color


########################### LOAD DATA ###########################
# (a) every bird carrying IBD with the founder: beak pointedness against the genome-wide fraction of SNPs in IBD
ibd_birds = load_ibd_beak_birds(load_master())

# (c) the scandens of the association test with at least half of the chromosome 6 region in IBD, or none, tested with the same LMM as the
# windows (chromosome 6 LOCO GRM, covariates sex, IBD outside chromosome 6, and PC(Body Size))
chrom_code = {name: code for code, name in load_chromosomes().items()}[CHROM]
birds, eig = region_test_birds(SPP, chrom_code, START, END)
beta, p = region_lmm(birds, eig, TRAIT)
print(f"chr{CHROM}:{START}-{END} {SPP} {TRAIT}: {birds['ibd_group'].value_counts().to_dict()}, beta = {beta:.6f}, p = {p:.6f}")

# (b) windows of the association test
genome = genome_layout()
results = load_results(genome)
manhattan = results[(results["species"] == SPP) & (results["trait"] == TRAIT)]


########################### CREATE PLOT ###########################
top_width = A_WIDTH + GAP_AC + C_WIDTH + GAP_CD + TREE_WIDTH + GROUPS_WIDTH + HEATMAP_WIDTH + LEGEND_WIDTH
b_width = top_width - 0.1  # ends under the heatmap's legends
b_height = b_width / B_ASPECT
width, height = LEFT + top_width + RIGHT, TOP + ROW_HEIGHT + ROW_GAP + b_height + BOTTOM
fig = plt.figure(figsize=(width, height))


def add_axes(left, bottom, w, h):
    """Axes placed in inches from the figure's bottom left corner."""
    return fig.add_axes([left / width, bottom / height, w / width, h / height])


row_bottom = BOTTOM + b_height + ROW_GAP
x = LEFT
ax_a = add_axes(x, row_bottom, A_WIDTH, ROW_HEIGHT)
x += A_WIDTH + GAP_AC
ax_c = add_axes(x, row_bottom, C_WIDTH, ROW_HEIGHT)
x += C_WIDTH + GAP_CD
ax_tree = add_axes(x, row_bottom, TREE_WIDTH, ROW_HEIGHT)
ax_groups = add_axes(x + TREE_WIDTH + 0.04, row_bottom, GROUPS_WIDTH, ROW_HEIGHT)
ax_d = add_axes(x + TREE_WIDTH + GROUPS_WIDTH + 0.08, row_bottom, HEATMAP_WIDTH, ROW_HEIGHT)
ax_b = add_axes(LEFT, BOTTOM, b_width, b_height)

plot_beak_vs_ibd(ax_a, ibd_birds, TRAIT, show_legend=True, stats_y=(1.075, 1.02))
style_legend(ax_a, loc="upper right")

plot_region_box(ax_c, birds, TRAIT, beta, p, SPP)
ax_c.set_title(f"{START / 1e6:g}-{END / 1e6:g} Mb on chr. {CHROM}", pad=30)
ax_c.texts[-1].set_text(ax_c.texts[-1].get_text().replace("Linear mixed model", "LMM"))  # the statistics under the title

# heatmap: the peak's position as a subtitle (as in B), and a compact legend (no islands, the scandens labels on two lines)
group_labels = {g: GROUPS.loc[g, "label"].replace(" (Daphne)", "").replace(" (Española)", "").replace(" (", "\n(") for g in GROUPS.index}
n_snps, n_shown, n_scandens = plot_peak_heatmap(ax_tree, ax_groups, ax_d, fontsize=HEATMAP_FONTSIZE, legend_y=(1, 0.62), peak_labels=False,
                                                group_labels=group_labels, colors=HEATMAP_COLORS)
peaks = pd.read_csv(PEAK_REGIONS, sep="\t")
ax_d.set_title(f"{n_snps} {HEATMAP_SNPS} SNPs in the FST peak{'s' if len(peaks) > 1 else ''}", pad=30)
ax_d.text(0.5, 1.02, ", ".join(f"{p['start'] / 1e6:.2f}-{p['end'] / 1e6:.2f}" for _, p in peaks.iterrows()) + f" Mb on chr. {CHROM}",
          transform=ax_d.transAxes, ha="center", va="bottom", fontsize="small", color="#333333")
print(f"heatmap: {n_shown} of the {n_scandens} scandens without IBD shown")

plot_manhattan(ax_b, manhattan, genome, SPP, "neg_log10_p", "$-log_{10}$p-value", show_threshold=True, point_size=8, tick_fontsize=13, combine_from="17")
ax_b.set_title(f"Association of {IBD_F} with {TRAIT} in " + species_label(SPP), loc="left", pad=16)

# gray band behind the chromosome 6 region in (b)
region_start, region_end = genome.loc[CHROM, "start"] + START, genome.loc[CHROM, "start"] + END
ax_b.axvspan(region_start, region_end, color="lightgray", alpha=0.5, linewidth=0, zorder=0)

# each letter over its panel's y axis label (D over the heatmap's tree), A and B (the bottom row) under each other
add_panel_letters(fig, [ax_a, ax_b, ax_c, ax_tree], align=[[0, 1]], anchor="label")
save_figure("main/figure4v2.png", fig, pdf=True)
