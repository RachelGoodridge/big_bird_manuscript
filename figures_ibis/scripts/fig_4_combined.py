import importlib.util
import sys
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ConnectionPatch
from common import ASSOCIATION, load_chromosomes, load_ibd_beak_birds, load_master, plot_beak_vs_ibd, save_figure, set_plot_theme, style_legend

# the IBD-phenotype association helpers (loaded by path since it is also named common.py)
_spec = importlib.util.spec_from_file_location("association_common", f"{ASSOCIATION}/scripts/common.py")
sys.modules["association_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["association_common"])
from association_common import CHR6_REGION, genome_layout, load_results, plot_manhattan, plot_region_box, region_lmm, region_test_birds

set_plot_theme()
SPP = "scandens"
TRAIT = "PC(Beak Pointedness)"
CHROM, START, END = CHR6_REGION

# where the arrow from the region in (c) ends, in (b)'s axes fractions: x from 0 (left edge of b's plot area) to 1 (right edge), y from 0 (bottom
# of the plot area, above the tick labels) to 1 (top); negative y is below the plot area (the tick labels end around -0.13)
ARROW_TIP = (0.3, -0.12)


########################### LOAD DATA ###########################
# (a) every bird carrying IBD with the founder: beak pointedness against the genome-wide fraction of SNPs in IBD
ibd_birds = load_ibd_beak_birds(load_master())

# (b) the scandens of the association test with at least half of the chromosome 6 region in IBD, or none, tested with the same LMM as the
# windows (chromosome 6 LOCO GRM, covariates sex, IBD outside chromosome 6, and PC(Body Size))
chrom_code = {name: code for code, name in load_chromosomes().items()}[CHROM]
birds, eig = region_test_birds(SPP, chrom_code, START, END)
beta, p = region_lmm(birds, eig, TRAIT)
print(f"chr{CHROM}:{START}-{END} {SPP} {TRAIT}: {birds['ibd_group'].value_counts().to_dict()}, beta = {beta:.6f}, p = {p:.6f}")

# (c) windows of the association test
genome = genome_layout()
results = load_results(genome)
manhattan = results[(results["species"] == SPP) & (results["trait"] == TRAIT)]


########################### CREATE PLOT ###########################
fig = plt.figure(figsize=(16, 12))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 0.75], hspace=0.32, wspace=0.25)
ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, :])

plot_beak_vs_ibd(ax_a, ibd_birds, TRAIT, show_legend=True, stats_y=(1.075, 1.02))
style_legend(ax_a, loc="upper right")

plot_region_box(ax_b, birds, TRAIT, beta, p, SPP)
ax_b.set_title(f"IBD on chr. {CHROM} ({START / 1e6:g}-{END / 1e6:g} Mb) in " + r"$\it{" + SPP + "}$", pad=30)

plot_manhattan(ax_c, manhattan, genome, SPP, "neg_log10_p", "$-log_{10}$p-value", show_threshold=True, point_size=8, tick_fontsize=11, combine_from="17")
ax_c.set_title(r"Association of $\mathrm{IBD_{F}}$ with " + f"{TRAIT} in " + r"$\it{" + SPP + "}$", loc="left", pad=16)  # left so the arrow from chr. 6 clears it

# gray band behind the chromosome 6 region in (c), with an arrow from it up to (b)
region_start, region_end = genome.loc[CHROM, "start"] + START, genome.loc[CHROM, "start"] + END
ax_c.axvspan(region_start, region_end, color="lightgray", alpha=0.5, linewidth=0, zorder=0)
band_top = mtransforms.blended_transform_factory(ax_c.transData, ax_c.transAxes)
fig.add_artist(ConnectionPatch(xyA=((region_start + region_end) / 2, 1.0), coordsA=band_top, xyB=ARROW_TIP, coordsB=ax_b.transAxes,
                               arrowstyle="-|>", color="darkgray", linewidth=1.5, mutation_scale=18))

save_figure("main/figure4.png", fig)
