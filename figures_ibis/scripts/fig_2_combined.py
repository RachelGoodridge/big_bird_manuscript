import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec
from common import (FOCAL_SPECIES, GROUPS, add_bill_pcs, add_freq_colorbar, add_panel_letters, chromosome_lengths, format_year_axis, founder_ibd_per_bird,
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
# (a) beak pointedness PC fit on every measured fortis / scandens on Daphne, yearly mean +/- sd (years with at least 15 birds of a species),
# with every point the same size as the yearly means in panels c-d (not scaled by the number of birds)
master = load_master()
measured = add_bill_pcs(load_measured_fortis_scandens(master))
mean, sd, _, n_scaled = yearly_trait_summary(measured, ["PC(Beak Pointedness)"])
point_size = pd.Series(25, index=n_scaled.index)

# (b-d) fraction of SNPs in IBD with the founder for every bird, and per year averages / carrier frequencies
birds = founder_ibd_per_bird(master)
summary = order_groups(yearly_summary(birds), FOCAL_SPECIES, by=["Year", "Species"])
birds = order_groups(birds, GROUPS)

# (e) IBD frequency along the genome in each year, on a color scale shared by both years and species
chr_dict = load_chromosomes()
chrom_lengths = chromosome_lengths(load_bim())
profiles = {year: ibd_frequency_profiles(master, year, chrom_lengths, BIN_SIZE) for year in YEARS}
vmax = max(fq.max() for p in profiles.values() for _, fq in p.values())


########################### CREATE PLOT ###########################
# ---------------- layout settings you may want to change (fractions of the figure, except the ratios) ----------------
LEFT, RIGHT, TOP, BOTTOM = 0.06, 0.97, 0.96, 0.05   # outer edges of the plot areas of panels a-d
COLUMN_GAP, ROW_GAP = 0.24, 0.30                    # blank space between panels, relative to the average panel width / height
ROW_HEIGHTS = [1, 1, 0.02, 0.45]                    # panels a-b, panels c-d, the legend, panel e
E_RAISE = 0.022                                     # panel e moved up by this much, closer to the legend
COLORBAR_GAP, COLORBAR_WIDTH = 0.012, 0.010          # between panel e and its colorbar, and the colorbar's width
# ----------------------------------------------------------------------------------------------------------------------
# 2x2 grid of panels a-d, then the legend of panels a-d in a row of its own, then the genome-wide IBD tracks (e) spanning the full width
fig = plt.figure(figsize=(20, 21))
gs = GridSpec(4, 2, figure=fig, height_ratios=ROW_HEIGHTS, left=LEFT, right=RIGHT, top=TOP, bottom=BOTTOM, wspace=COLUMN_GAP, hspace=ROW_GAP)
ax_a, ax_b, ax_c, ax_d = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]
ax_e = fig.add_subplot(gs[3, :])
box_e = ax_e.get_position()
ax_e.set_position([box_e.x0, box_e.y0 + E_RAISE, box_e.width, box_e.height])  # its width is set below, once its labels can be measured

plot_trait_means(ax_a, mean, sd, point_size, "PC(Beak Pointedness)", POINT_SCALE, LINE_SCALE)
ax_a.set_ylabel("Mean PC(Beak Pointedness)")
format_year_axis(ax_a)
ax_a.tick_params(axis="x", rotation=30)
sns.despine(top=True, right=True, left=False, bottom=False, ax=ax_a)

handles, labels = plot_ibd_over_time([ax_b, ax_c, ax_d], birds, summary, POINT_SCALE, LINE_SCALE)

plot_ibd_track(fig, ax_e, profiles, YEARS, chr_dict, chrom_lengths, vmax, BIN_SIZE, fontsize=17, combine_from="17")

# one legend for panels a-d in a row under them, so it covers no data (panel b has every group)
box_legend = gs[2, :].get_position(fig)
leg = fig.legend(handles, labels, loc="center", ncol=len(labels), bbox_to_anchor=(0.5, (box_legend.y0 + box_legend.y1) / 2), frameon=False, fontsize=22,
                 handletextpad=0.1)
italicize_species_legend(leg)

# panel e, with its labels on the left and its colorbar and the colorbar's label on the right, is made exactly as wide as panels a-b with
# theirs: its text starts where panel a's y axis label starts, and the colorbar's label ends where panel b's plot area ends. Everything is
# drawn once to measure how far the labels and the colorbar stick out of panel e, then panel e is narrowed by that much
def colorbar_box(box):
    return [box.x1 + COLORBAR_GAP, box.y0, COLORBAR_WIDTH, box.height]


def figure_box(artist):
    return artist.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())


bbox_e = ax_e.get_position()
colorbar = add_freq_colorbar(fig, colorbar_box(bbox_e), "Reds", vmax, fontsize=18)
fig.canvas.draw()
left_overhang = bbox_e.x0 - figure_box(ax_e).x0
right_overhang = figure_box(colorbar.ax).x1 - bbox_e.x1
left_edge, right_edge = figure_box(ax_a).x0, ax_b.get_position().x1
ax_e.set_position([left_edge + left_overhang, bbox_e.y0, right_edge - right_overhang - left_edge - left_overhang, bbox_e.height])
colorbar.ax.set_position(colorbar_box(ax_e.get_position()))

# the letters start over the panels' y axis labels, lined up within each column (A, C, E and B, D); larger, like the rest of this figure's
# text
add_panel_letters(fig, [ax_a, ax_b, ax_c, ax_d, ax_e], fontsize=32, pad=16, align=[[0, 2, 4], [1, 3]], anchor="label")
save_figure("main/figure2.png", fig, pdf=True)
