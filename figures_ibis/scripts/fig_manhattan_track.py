import matplotlib.pyplot as plt
from common import (add_freq_colorbar, chromosome_lengths, ibd_frequency_profiles, load_bim, load_chromosomes, load_master, plot_ibd_track,
                    save_figure, set_genome_plot_theme)

set_genome_plot_theme()
YEARS = [2002, 2022]  # top to bottom within each species
BIN_SIZE = 0.1


########################### LOAD DATA ###########################
chr_dict = load_chromosomes()
chrom_lengths = chromosome_lengths(load_bim())
master = load_master()

# IBD frequency along the genome in each year, on a color scale shared by both years and species
profiles = {year: ibd_frequency_profiles(master, year, chrom_lengths, BIN_SIZE) for year in YEARS}
vmax = max(fq.max() for p in profiles.values() for _, fq in p.values())
print(f"Max IBD frequency (color scale): {vmax:.4f}")


########################### MANHATTAN-STYLE FREQUENCY TRACK PLOT ###########################
# 4 stacked genome-wide tracks (top→bottom: fortis 2002, fortis 2022, scandens 2002, scandens 2022) colored by IBD segment frequency
fig, ax = plt.subplots(figsize=(30, 4))
plot_ibd_track(fig, ax, profiles, YEARS, chr_dict, chrom_lengths, vmax, BIN_SIZE)

# a single colorbar shared by both species, just right of the tracks
fig.canvas.draw()
bbox_main = ax.get_position()
add_freq_colorbar(fig, [bbox_main.x1 + 0.010, bbox_main.y0, 0.010, bbox_main.height], "Reds", vmax)

save_figure("main/IBD_manhattan_track_2022_2002.png", fig)
