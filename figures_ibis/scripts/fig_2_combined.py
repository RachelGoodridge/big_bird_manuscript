import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from common import (FOCAL_SPECIES, GROUPS, add_bill_pcs, add_freq_colorbar, chromosome_lengths, format_year_axis, founder_ibd_per_bird,
                    ibd_frequency_profiles, italicize_species_legend, load_bim, load_chromosomes, load_master, load_measured_fortis_scandens,
                    order_groups, plot_ibd_over_time, plot_ibd_track, plot_trait_means, save_figure, set_plot_theme, yearly_summary,
                    yearly_trait_summary)

set_plot_theme()
# larger text than the shared theme, since this figure is shrunk more to fit the page (only affects this script)
plt.rcParams.update({"axes.labelsize": 22, "xtick.labelsize": 18, "ytick.labelsize": 18})
POINT_SCALE = 3  # larger points and thicker lines / error bars than the shared figures
LINE_SCALE = 2
YEARS = [2002, 2022]  # IBD tracks, top to bottom within each species
BIN_SIZE = 0.1


########################### LOAD DATA ###########################
# (a-c) fraction of SNPs in IBD with the founder for every bird, and per year averages / carrier frequencies
master = load_master()
birds = founder_ibd_per_bird(master)
summary = order_groups(yearly_summary(birds), FOCAL_SPECIES, by=["Year", "Species"])
birds = order_groups(birds, GROUPS)

# (d) beak pointedness PC fit on every measured fortis / scandens on Daphne, yearly mean +/- sd (years with at least 15 birds of a species)
measured = add_bill_pcs(load_measured_fortis_scandens(master))
mean, sd, _, point_size = yearly_trait_summary(measured, ["PC(Beak Pointedness)"])

# (e) IBD frequency along the genome in each year, on a color scale shared by both years and species
chr_dict = load_chromosomes()
chrom_lengths = chromosome_lengths(load_bim())
profiles = {year: ibd_frequency_profiles(master, year, chrom_lengths, BIN_SIZE) for year in YEARS}
vmax = max(fq.max() for p in profiles.values() for _, fq in p.values())


########################### CREATE PLOT ###########################
# 2x2 grid of panels a-d, with the genome-wide IBD tracks (e) spanning the full width underneath
fig = plt.figure(figsize=(20, 20))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 0.45])
ax_a, ax_b, ax_c, ax_d = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]
ax_e = fig.add_subplot(gs[2, :])

handles, labels = plot_ibd_over_time([ax_a, ax_b, ax_c], birds, summary, POINT_SCALE, LINE_SCALE)

plot_trait_means(ax_d, mean, sd, point_size, "PC(Beak Pointedness)", POINT_SCALE, LINE_SCALE)
ax_d.set_ylabel("Mean PC(Beak Pointedness)")
format_year_axis(ax_d)
ax_d.tick_params(axis="x", rotation=30)
sns.despine(top=True, right=True, left=False, bottom=False, ax=ax_d)

plot_ibd_track(fig, ax_e, profiles, YEARS, chr_dict, chrom_lengths, vmax, BIN_SIZE, fontsize=17, combine_from="17")

# one legend for panels a-d in a row above them, so it covers no data (panel a has every group)
leg = fig.legend(handles, labels, loc="upper center", ncol=len(labels), bbox_to_anchor=(0.5, 0.995), frameon=False, fontsize=22)
italicize_species_legend(leg)

# leave room on the right for the IBD frequency colorbar of panel e
plt.tight_layout(rect=(0, 0, 0.94, 0.97))
bbox_e = ax_e.get_position()
add_freq_colorbar(fig, [bbox_e.x1 + 0.012, bbox_e.y0, 0.010, bbox_e.height], "Reds", vmax, fontsize=18)

save_figure("main/figure2.png", fig)
