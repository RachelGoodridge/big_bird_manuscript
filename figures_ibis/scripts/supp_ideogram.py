import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec
from common import (FOCAL_SPECIES, add_freq_colorbar, chromosome_lengths, ibd_frequency_profiles,
                    load_bim, load_chromosomes, load_master, save_figure, set_genome_plot_theme)

set_genome_plot_theme()
year_choice = 2022
year_compare = 2002  # only used so the color scale matches the Manhattan track figure
BIN_SIZE = 0.1


########################### LOAD DATA ###########################
chr_dict = load_chromosomes()
chroms_ordered = list(chr_dict.keys())
chrom_lengths = chromosome_lengths(load_bim())
master = load_master()

# IBD frequency along the genome, on a color scale shared by both years and species
profiles = ibd_frequency_profiles(master, year_choice, chrom_lengths, BIN_SIZE)
vmax = max(fq.max() for p in [profiles, ibd_frequency_profiles(master, year_compare, chrom_lengths, BIN_SIZE)] for _, fq in p.values())
print(f"Max IBD frequency (color scale): {vmax:.4f}")


########################### DRAWING HELPER ###########################
PLOT_HEIGHT = 3.0
CHROM_WIDTH = 0.55
X_SPACING   = 1.3
# cap_h target = CHROM_WIDTH * CORNER_FRAC. Large chromosomes get a flat semi-ellipse cap;
# small chromosomes are clamped to chrom_height/2, giving a full pill/ellipse shape.
CORNER_FRAC = 0.125

def make_stadium_path(x0, y0, w, h, cap_h):
    """
    Return a Path for a stadium shape: straight vertical sides joined to semi-elliptical
    top and bottom caps. The ellipse tangent is vertical at each side-cap junction so
    the shape is smooth with no sharp corners. cap_h is the cap's vertical semi-axis,
    clamped to h/2 so the two caps meet cleanly on very short chromosomes.
    """
    cap_h = min(cap_h, h / 2)
    a  = w / 2            # horizontal semi-axis of the elliptical cap
    k  = 0.5523           # bezier constant for a quarter-ellipse approximation
    xc = x0 + a           # horizontal centre
    bot_y = y0 + cap_h
    top_y = y0 + h - cap_h
    path_data = [
        (Path.MOVETO,    (x0,          bot_y)),
        # bottom semi-ellipse: left side → nadir → right side (2 quarter-ellipses)
        (Path.CURVE4,    (x0,          bot_y - k * cap_h)),
        (Path.CURVE4,    (xc - k * a,  y0)),
        (Path.CURVE4,    (xc,          y0)),
        (Path.CURVE4,    (xc + k * a,  y0)),
        (Path.CURVE4,    (x0 + w,      bot_y - k * cap_h)),
        (Path.CURVE4,    (x0 + w,      bot_y)),
        # right straight side
        (Path.LINETO,    (x0 + w,      top_y)),
        # top semi-ellipse: right side → apex → left side
        (Path.CURVE4,    (x0 + w,      top_y + k * cap_h)),
        (Path.CURVE4,    (xc + k * a,  y0 + h)),
        (Path.CURVE4,    (xc,          y0 + h)),
        (Path.CURVE4,    (xc - k * a,  y0 + h)),
        (Path.CURVE4,    (x0,          top_y + k * cap_h)),
        (Path.CURVE4,    (x0,          top_y)),
        # left straight side back to start
        (Path.LINETO,    (x0,          bot_y)),
        (Path.CLOSEPOLY, (x0,          bot_y)),
    ]
    codes, verts = zip(*path_data)
    return Path(verts, codes)


def draw_chromosome(ax, chrom_len_mb, max_len_mb, bin_centers, freq, x_center, highlight_cmap, sp_norm, bin_size=BIN_SIZE):
    chrom_height = (chrom_len_mb / max_len_mb) * PLOT_HEIGHT
    x0 = x_center - CHROM_WIDTH / 2
    # Clamping at h/9 guarantees ≥50% straight sides for every chromosome size,
    # keeping the stadium shape consistent across the top and bottom rows.
    cap_h = min(chrom_height / 9.0, CHROM_WIDTH * CORNER_FRAC)
    body_path = make_stadium_path(x0, 0.0, CHROM_WIDTH, chrom_height, cap_h)

    # grey body
    ax.add_patch(mpatches.PathPatch(body_path, facecolor="#d8d8d8", edgecolor="#555555", linewidth=0.8, zorder=2))

    # coloured bands clipped to the stadium shape
    if freq.max() > 0:
        cmap = plt.get_cmap(highlight_cmap)
        # Build a clip patch from the same path (must be a separate patch instance)
        clip_patch = mpatches.PathPatch(body_path, transform=ax.transData)
        ax.add_patch(clip_patch)
        clip_patch.set_visible(False)

        for bc, fq in zip(bin_centers, freq):
            if fq <= 0:
                continue
            y_bot = max((bc - bin_size / 2) / chrom_len_mb * chrom_height, 0.0)
            y_top = min((bc + bin_size / 2) / chrom_len_mb * chrom_height, chrom_height)
            band_h = y_top - y_bot
            if band_h <= 0:
                continue
            band = Rectangle((x0, y_bot), CHROM_WIDTH, band_h, linewidth=0, facecolor=cmap(sp_norm(fq)), zorder=3)
            band.set_clip_path(clip_patch)
            ax.add_patch(band)

    # outline on top (same stadium path)
    ax.add_patch(mpatches.PathPatch(make_stadium_path(x0, 0.0, CHROM_WIDTH, chrom_height, cap_h), facecolor="none", edgecolor="#555555", linewidth=0.8, zorder=4))
    return chrom_height


########################### FIGURE LAYOUT ###########################
SPECIES_CONFIG = {
    "fortis":   {"cmap": "Blues",  "title": r"$\it{G.\ fortis}$ " + f"({year_choice})"},
    "scandens": {"cmap": "Greens", "title": r"$\it{G.\ scandens}$ " + f"({year_choice})"},
}

top_chroms    = chroms_ordered[:15]
bottom_chroms = chroms_ordered[15:30]
max_len_mb    = chrom_lengths.max()
top_max_mb    = max(chrom_lengths[c] for c in top_chroms)
bottom_max_mb = max(chrom_lengths[c] for c in bottom_chroms)

fig = plt.figure(figsize=(30, 13))

# GridSpec: 2 species columns, 2 chromosome rows.
# Colorbars placed manually after draw via bbox so they never overlap.
# left/right margins leave room for colorbars outside both species columns.
gs = GridSpec(
    2, 2,
    figure=fig,
    height_ratios=[top_max_mb, bottom_max_mb * 6.5],  # bottom row much taller
    width_ratios=[1, 1],
    hspace=0.00,    # minimal vertical gap between rows
    wspace=0.30,    # wider gap so fortis cbar+label fits before scandens panel starts
    top=0.91,
    bottom=0.04,
    left=0.05,
    right=0.95,
)
axes_grid = {(row, species): fig.add_subplot(gs[row, col]) for col, species in enumerate(FOCAL_SPECIES) for row in [0, 1]}

FS_CHROM_LABEL = 13    # chromosome name below pill
FS_YTICK       = 13    # y-axis tick numbers
FS_YLABEL      = 14    # "Chromosome length (Mb)"
FS_TITLE       = 17    # species subtitle (G. fortis / G. scandens)

sp_norm = Normalize(vmin=0, vmax=vmax)
for species, cfg in SPECIES_CONFIG.items():
    for row_idx, chrom_group in enumerate([top_chroms, bottom_chroms]):
        ax = axes_grid[(row_idx, species)]
        row_max_mb = max(chrom_lengths[c] for c in chrom_group)
        ax.set_xlim(-CHROM_WIDTH * 0.8, len(chrom_group) * X_SPACING + CHROM_WIDTH * 0.3)

        for i, chrom in enumerate(chrom_group):
            bc, fq = profiles[(species, chrom)]
            x_center = i * X_SPACING + CHROM_WIDTH / 2
            draw_chromosome(ax, chrom_lengths[chrom], max_len_mb, bc, fq, x_center, cfg["cmap"], sp_norm)
            ax.text(x_center, -0.13, str(chr_dict[chrom]), ha="center", va="top", fontsize=FS_CHROM_LABEL, color="#333333")

        # y-axis: spine and ticks clipped to data range
        if row_idx == 1:
            # Bottom row: round up to next multiple of 5, step every 5 Mb
            y_max_tick_mb  = int(np.ceil(row_max_mb / 5)) * 5
            y_tick_mb_vals = np.arange(0, y_max_tick_mb + 1, 5)
        else:
            y_tick_mb_vals = np.arange(0, int(row_max_mb) + 1, max(1, int(row_max_mb // 5)))
        y_tick_positions = [(v / max_len_mb) * PLOT_HEIGHT for v in y_tick_mb_vals]
        # ylim must encompass both the tallest chromosome and the last tick mark
        y_top_data = max((row_max_mb / max_len_mb) * PLOT_HEIGHT, y_tick_positions[-1])
        ax.set_ylim(-0.45, y_top_data + 0.15)
        ax.set_yticks(y_tick_positions)
        ax.set_yticklabels([f"{int(v)}" for v in y_tick_mb_vals], fontsize=FS_YTICK)
        ax.yaxis.set_visible(True)
        ax.spines[["top", "right", "bottom"]].set_visible(False)
        ax.spines["left"].set_visible(True)
        ax.spines["left"].set_color("#aaaaaa")
        ax.spines["left"].set_bounds(y_tick_positions[0], y_tick_positions[-1])
        ax.tick_params(axis="y", length=3, color="#aaaaaa")
        ax.xaxis.set_visible(False)

        # y-label on ALL four subplots, centered on tick midpoint not full ylim
        y_label_pos = (y_tick_positions[0] + y_tick_positions[-1]) / 2
        ax.set_ylabel("Chromosome length (Mb)", fontsize=FS_YLABEL, labelpad=8)
        ax.yaxis.set_label_coords(-0.04, (y_label_pos - ax.get_ylim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0]))

        if row_idx == 0:
            ax.set_title(cfg["title"], fontsize=FS_TITLE, fontfamily="DejaVu Serif", pad=6, color="#111111")

# colorbars just right of each species' top-row panel, matching its height
fig.canvas.draw()
for species, cfg in SPECIES_CONFIG.items():
    bbox_top = axes_grid[(0, species)].get_position()
    add_freq_colorbar(fig, [bbox_top.x1 + 0.008, bbox_top.y0, 0.010, bbox_top.height], cfg["cmap"], vmax)

save_figure(f"supplement/IBD_ideogram_{year_choice}.png")
