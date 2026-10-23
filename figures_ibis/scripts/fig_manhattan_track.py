import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
import matplotlib.transforms as mtransforms
from common import add_freq_colorbar, chromosome_lengths, ibd_frequency_profiles, load_bim, load_chromosomes, load_master, save_figure, set_genome_plot_theme

set_genome_plot_theme()
year_choice = 2022
year_choice2 = 2002
BIN_SIZE = 0.1


########################### LOAD DATA ###########################
chr_dict = load_chromosomes()
chroms_ordered = list(chr_dict.keys())
chrom_lengths = chromosome_lengths(load_bim())
master = load_master()

# IBD frequency along the genome in each year, on a color scale shared by both years and species
profiles = {year: ibd_frequency_profiles(master, year, chrom_lengths, BIN_SIZE) for year in [year_choice, year_choice2]}
vmax = max(fq.max() for p in profiles.values() for _, fq in p.values())
print(f"Max IBD frequency (color scale): {vmax:.4f}")


########################### MANHATTAN-STYLE FREQUENCY TRACK PLOT ###########################
# Chromosomes laid end-to-end (Manhattan-plot x-axis); 4 stacked genome-wide tracks
# (top→bottom: fortis year_choice2, fortis year_choice, scandens year_choice2, scandens year_choice)
# colored by IBD segment frequency. Gap scales with the largest chromosome so tiny
# chromosomes stay visually distinct without wasting space on huge ones. Wide enough that
# the white axes background reads as a clear boundary instead of a thin separator line.
CHROM_GAP_MB = chrom_lengths.max() * 0.02

chrom_offset = {}
cum_pos = 0.0
for chrom in chroms_ordered:
    chrom_offset[chrom] = cum_pos
    cum_pos += chrom_lengths[chrom] + CHROM_GAP_MB
total_genome_mb = cum_pos - CHROM_GAP_MB

TRACK_HEIGHT = 0.5
ROW_GAP      = 0.1   # vertical gap between rows, in TRACK_HEIGHT units
ROW_PITCH    = TRACK_HEIGHT + ROW_GAP
TRACK_ROWS   = [("fortis", year_choice2), ("fortis", year_choice), ("scandens", year_choice2), ("scandens", year_choice)]
n_rows = len(TRACK_ROWS)
FS_LABEL = 13

fig, ax = plt.subplots(figsize=(30, 4))
cmap = plt.get_cmap("Reds")
norm = Normalize(vmin=0, vmax=vmax)

row_y0 = {}
for row_idx, (species, yr) in enumerate(TRACK_ROWS):
    y0 = (n_rows - 1 - row_idx) * ROW_PITCH
    row_y0[row_idx] = y0

    for chrom in chroms_ordered:
        x_off = chrom_offset[chrom]
        ax.add_patch(Rectangle((x_off, y0), chrom_lengths[chrom], TRACK_HEIGHT, facecolor="#d8d8d8", edgecolor="none", zorder=1))
        bc, fq = profiles[yr][(species, chrom)]
        for bcenter, f in zip(bc, fq):
            if f <= 0:
                continue
            x_lo  = x_off + max(bcenter - BIN_SIZE / 2, 0.0)
            x_hi  = x_off + min(bcenter + BIN_SIZE / 2, chrom_lengths[chrom])
            width = x_hi - x_lo
            if width <= 0:
                continue
            ax.add_patch(Rectangle((x_lo, y0), width, TRACK_HEIGHT, facecolor=cmap(norm(f)), edgecolor="none", zorder=2))

# chromosome and row boundaries are conveyed by the gaps themselves (white axes
# background shows through between chromosome rectangles and between row rectangles)

# thin outline box around each species' pair of rows (2002 + 2022 tracks)
species_row_indices = {}
for row_idx, (species, yr) in enumerate(TRACK_ROWS):
    species_row_indices.setdefault(species, []).append(row_idx)

for species, rows in species_row_indices.items():
    block_y0 = row_y0[max(rows)]
    block_y1 = row_y0[min(rows)] + TRACK_HEIGHT
    ax.add_patch(Rectangle((0, block_y0), total_genome_mb, block_y1 - block_y0, facecolor="none", edgecolor="#333333", linewidth=1.0, zorder=5, clip_on=False))

ax.set_xlim(0, total_genome_mb)
ax.set_ylim(0, (n_rows - 1) * ROW_PITCH + TRACK_HEIGHT)

# individual labels through chromosome "19"; the tail of tiny micro-chromosomes
# (21, 22, ..., 28, LGE22) is too narrow for separate labels, so it collapses
# into one combined tick, mirroring the grouped-tail style of the reference Manhattan plot
split_idx = [chr_dict[c] for c in chroms_ordered].index("21")
individual_chroms = chroms_ordered[:split_idx]
combined_chroms   = chroms_ordered[split_idx:]
combined_x0       = chrom_offset[combined_chroms[0]]
combined_x1       = chrom_offset[combined_chroms[-1]] + chrom_lengths[combined_chroms[-1]]

x_tick_positions = [chrom_offset[c] + chrom_lengths[c] / 2 for c in individual_chroms]
x_tick_labels    = [str(chr_dict[c]) for c in individual_chroms]
x_tick_positions.append((combined_x0 + combined_x1) / 2)
x_tick_labels.append("21-28, LGE22")
ax.set_xticks(x_tick_positions)
ax.set_xticklabels(x_tick_labels, fontsize=FS_LABEL)
ax.set_xlabel("Chromosome", fontsize=14, labelpad=8)

# unlabeled minor ticks: one centered in each inter-chromosome gap, plus one at
# the very start and one at the very end of the whole genome
boundary_ticks = [0.0]
for chrom_a, chrom_b in zip(chroms_ordered[:-1], chroms_ordered[1:]):
    gap_start = chrom_offset[chrom_a] + chrom_lengths[chrom_a]
    gap_end   = chrom_offset[chrom_b]
    boundary_ticks.append((gap_start + gap_end) / 2)
boundary_ticks.append(total_genome_mb)
ax.set_xticks(boundary_ticks, minor=True)
ax.tick_params(axis="x", which="minor", length=3, color="#aaaaaa", width=0.8)
ax.tick_params(axis="x", which="major", length=0)

row_y_centers = [(n_rows - 1 - i) * ROW_PITCH + TRACK_HEIGHT / 2 for i in range(n_rows)]
ax.set_yticks([])

# two label columns hugging the left edge of the axis, both rotated vertical
# like the colorbar label: species name spanning its two rows in the outer
# column, year (per row) in the inner column closer to the axis. Offsets are
# in points (not axes fraction) so the gap stays small and fixed regardless
# of the axes' width.
SPECIES_OFFSET_PT = 30
YEAR_OFFSET_PT     = 10

base_trans    = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
species_trans = mtransforms.offset_copy(base_trans, fig=fig, x=-SPECIES_OFFSET_PT, y=0, units="points")
year_trans    = mtransforms.offset_copy(base_trans, fig=fig, x=-YEAR_OFFSET_PT, y=0, units="points")

for species, rows in species_row_indices.items():
    y_center = (row_y_centers[min(rows)] + row_y_centers[max(rows)]) / 2
    ax.text(0.0, y_center, r"$\it{" + species + "}$", transform=species_trans, ha="center", va="center", rotation=90, fontsize=FS_LABEL, clip_on=False)

for row_idx, (species, yr) in enumerate(TRACK_ROWS):
    ax.text(0.0, row_y_centers[row_idx], f"{yr}", transform=year_trans, ha="center", va="center", rotation=90, fontsize=FS_LABEL, clip_on=False)

ax.spines[["top", "right", "left"]].set_visible(False)
ax.spines["bottom"].set_color("#aaaaaa")
ax.spines["bottom"].set_linewidth(0.8)
ax.tick_params(axis="both", length=0)

# a single colorbar shared by both species, just right of the tracks
fig.canvas.draw()
bbox_main = ax.get_position()
add_freq_colorbar(fig, [bbox_main.x1 + 0.010, bbox_main.y0, 0.010, bbox_main.height], "Reds", vmax)

save_figure(f"main/IBD_manhattan_track_{year_choice}_{year_choice2}.png", fig)
