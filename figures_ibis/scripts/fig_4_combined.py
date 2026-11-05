import importlib.util
import sys
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import seaborn as sns
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ConnectionPatch
from common import ASSOCIATION, COLOR_MAP, load_chromosomes, load_ibd_beak_birds, load_master, plot_beak_vs_ibd, save_figure, set_plot_theme, style_legend

# the IBD-phenotype association helpers (loaded by path since it is also named common.py)
_spec = importlib.util.spec_from_file_location("association_common", f"{ASSOCIATION}/scripts/common.py")
sys.modules["association_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["association_common"])
from association_common import chromosome_lmm, chromosome_test_birds, genome_layout, load_results, plot_manhattan

set_plot_theme()
SPP = "scandens"
TRAIT = "PC(Beak Pointedness)"
CHROM = "6"


########################### LOAD DATA ###########################
# (a) every bird carrying IBD with the founder: beak pointedness against the genome-wide fraction of SNPs in IBD
ibd_birds = load_ibd_beak_birds(load_master())

# (b) the scandens of the association test: beak pointedness against the fraction of chromosome 6 SNPs in IBD, tested with the same LMM as
# the windows (chromosome 6 LOCO GRM, covariates sex, IBD outside chromosome 6, and PC(Body Size)), but with continuous IBD over the whole chromosome
chrom_code = {name: code for code, name in load_chromosomes().items()}[CHROM]
birds, eig = chromosome_test_birds(SPP, chrom_code)
beta, p = chromosome_lmm(birds, eig, TRAIT)
carriers, non_carriers = birds[birds["chr_ibd"] > 0], birds[birds["chr_ibd"] == 0]
print(f"chr{CHROM} {SPP} {TRAIT}: n = {len(birds)}, {(birds['chr_ibd'] > 0).sum()} with IBD on chr{CHROM}, beta = {beta:.6f}, p = {p:.6f}")

# (c) windows of the association test
genome = genome_layout()
results = load_results(genome)
manhattan = results[(results["species"] == SPP) & (results["trait"] == TRAIT)]


########################### CREATE PLOT ###########################
fig = plt.figure(figsize=(16, 12))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1, 0.75], hspace=0.45, wspace=0.25)
ax_a, ax_b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, :])

plot_beak_vs_ibd(ax_a, ibd_birds, TRAIT, show_legend=True, stats_y=(1.075, 1.02))
style_legend(ax_a, loc="upper right")

# (b) the birds carrying IBD on chromosome 6 (as in a, birds without IBD are not shown), with the mean of the birds without any as a reference
sns.scatterplot(x=carriers["chr_ibd"], y=carriers[TRAIT], color=COLOR_MAP[SPP], s=25, alpha=0.75, ax=ax_b, zorder=2)
ax_b.axhline(non_carriers[TRAIT].mean(), color="darkgray", linestyle="--", linewidth=1.5, zorder=1)
ax_b.text(0.17, non_carriers[TRAIT].mean() - 0.03, f"mean without IBD (n = {len(non_carriers)})", ha="left", va="top", fontsize="x-small",
          color="dimgray")
ax_b.text(0.03, 1.02, rf"$\it{{{SPP}}}$ GRM-LMM: β = {beta:.3f}, p = {p:.3f}", transform=ax_b.transAxes, va="bottom", ha="left", fontsize="small",
          color=COLOR_MAP[SPP])
ax_b.set_xlabel(f"Proportion of Chr. {CHROM} SNPs in IBD")
ax_b.set_ylabel(TRAIT)
ax_b.set_xlim(-0.02, 1.02)
sns.despine(top=True, right=True, left=True, bottom=True, ax=ax_b)

plot_manhattan(ax_c, manhattan, genome, SPP, "neg_log10_p", "$-log_{10}$p-value", show_threshold=True, point_size=8, tick_fontsize=11, combine_from="17")
ax_c.set_title(f"{TRAIT} in " + r"$\it{" + SPP + "}$", pad=30)

# bracket over chromosome 6 in (c), with an arrow up to (b)
start, length = genome.loc[CHROM, "start"], genome.loc[CHROM, "length"]
bracket = mtransforms.blended_transform_factory(ax_c.transData, ax_c.transAxes)
y0, y1 = 1.02, 1.05
ax_c.plot([start, start, start + length, start + length], [y0, y1, y1, y0], transform=bracket, color="darkgray", linewidth=1.5, clip_on=False)
fig.add_artist(ConnectionPatch(xyA=(start + length / 2, y1), coordsA=bracket, xyB=(0.38, -0.125), coordsB=ax_b.transAxes,
                               arrowstyle="-|>", color="darkgray", linewidth=1.5, mutation_scale=18))

save_figure("main/figure4.png", fig)
