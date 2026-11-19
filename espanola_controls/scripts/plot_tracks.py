import importlib.util
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from common import CONTROLS, FOUNDER, MERGE_GAP_CM, MIN_LENGTH, ORIGINAL, REPO, load_focal_birds, load_master, original_min_length

# where along the genome the founder and each Española conirostris share IBD with Daphne birds, like figure 2E: one row per focal bird with all
# years combined, colored by the fraction of Daphne fortis / scandens sharing an IBD segment with it in each 0.1 Mb bin. The segments are the
# filtered, trimmed ones (output/trim_seg_espanola.pkl; with --min-length, those of the minimum length test)
_spec = importlib.util.spec_from_file_location("figures_common", f"{REPO}/figures_ibis/scripts/common.py")
sys.modules["figures_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["figures_common"])
from figures_common import FOCAL_SPECIES, add_freq_colorbar, chromosome_lengths, load_bim, load_chromosomes, set_genome_plot_theme, species_label  # noqa: E402

TEST = "--min-length" in sys.argv[1:]
OUT_DIR = MIN_LENGTH if TEST else CONTROLS
SEGMENTS = f"{MIN_LENGTH}/trim_seg_espanola.pkl" if TEST else f"{CONTROLS}/output/trim_seg_espanola.pkl"
FILTER_NOTE = f"segments > {original_min_length():.1f} cM (gaps ≤ {MERGE_GAP_CM:g} cM merged)" if TEST else ""
BIN_SIZE = 0.1  # Mb, as figure 2E

set_genome_plot_theme()
master = load_master().set_index("meaningful.unique")
chr_dict = load_chromosomes()
chroms = list(chr_dict.keys())
chrom_lengths = chromosome_lengths(load_bim())
focal_birds = load_focal_birds()

# Daphne fortis / scandens of the original run (the denominators), and the focal birds' segments with them (positions in Mb)
daphne = pd.Series(list(pd.read_pickle(f"{ORIGINAL}/output/coef_all_founder.pkl")["ind"]) + [FOUNDER])
species = daphne.map(master["Species"])
n_birds = species[species.isin(FOCAL_SPECIES)].value_counts()
seg = pd.read_pickle(SEGMENTS)
seg = seg[seg["partner"].isin(daphne)].copy()
seg["Species"] = seg["partner"].map(master["Species"])
seg = seg[seg["Species"].isin(FOCAL_SPECIES)]
seg["start_bin"] = (seg["phys_start_pos"] / BIN_SIZE).astype(int)  # positions in Mb
seg["end_bin"] = (seg["phys_end_pos"] / BIN_SIZE).astype(int)


def frequency(rows):
    """Fraction of a species' birds with a segment overlapping each bin, per chromosome (each bird counted once per bin)."""
    n_bins = rows["end_bin"] - rows["start_bin"] + 1
    bins = pd.DataFrame({"partner": np.repeat(rows["partner"].to_numpy(), n_bins), "chrom": np.repeat(rows["chrom"].to_numpy(), n_bins),
                         "bin": np.concatenate([np.arange(a, b + 1) for a, b in zip(rows["start_bin"], rows["end_bin"])]) if len(rows) else []})
    return bins.drop_duplicates().groupby(["chrom", "bin"]).size()


freqs = {(spp, focal): frequency(seg[(seg["focal"] == focal) & (seg["Species"] == spp)]) / n_birds[spp] for spp in FOCAL_SPECIES for focal in focal_birds}
vmax = max(f.max() for f in freqs.values() if len(f))
print(f"highest frequency {vmax:.3f}; per focal bird:")
print(pd.DataFrame({spp: [freqs[(spp, f)].max() if len(freqs[(spp, f)]) else 0 for f in focal_birds] for spp in FOCAL_SPECIES}, index=focal_birds).round(3).to_string())

# chromosomes end to end with gaps, as figure 2E
gap = chrom_lengths.max() * 0.02
offset = dict(zip(chroms, np.concatenate([[0], np.cumsum([chrom_lengths[c] + gap for c in chroms[:-1]])])))
total = offset[chroms[-1]] + chrom_lengths[chroms[-1]]

TRACK, PITCH, BLOCK_GAP = 0.5, 0.6, 0.5  # track height, row pitch, extra gap between the species blocks
rows = [(spp, focal) for spp in FOCAL_SPECIES for focal in focal_birds]
row_y = {}
y = 0.0
for spp in reversed(FOCAL_SPECIES):
    for focal in reversed(focal_birds):
        row_y[(spp, focal)] = y
        y += PITCH
    y += BLOCK_GAP

fig, ax = plt.subplots(figsize=(22, 13))
ax.set_rasterization_zorder(3)
cmap, norm = plt.get_cmap("Reds"), Normalize(vmin=0, vmax=vmax)
for (spp, focal), y0 in row_y.items():
    f = freqs[(spp, focal)]
    for chrom in chroms:
        ax.add_patch(Rectangle((offset[chrom], y0), chrom_lengths[chrom], TRACK, facecolor="#d8d8d8", edgecolor="none", zorder=1))
        if chrom not in f.index.get_level_values("chrom"):
            continue
        for b, value in f.loc[chrom].items():
            x0 = b * BIN_SIZE
            if x0 >= chrom_lengths[chrom]:
                continue
            ax.add_patch(Rectangle((offset[chrom] + x0, y0), min(BIN_SIZE, chrom_lengths[chrom] - x0), TRACK, facecolor=cmap(norm(value)),
                                   edgecolor="none", zorder=2))

# outline each species block
for spp in FOCAL_SPECIES:
    ys = [row_y[(spp, f)] for f in focal_birds]
    ax.add_patch(Rectangle((0, min(ys)), total, max(ys) + TRACK - min(ys), facecolor="none", edgecolor="#333333", linewidth=1.0, zorder=5))

# labels: each focal bird left of its row, the species left of its block
FONT = 15
row_trans = mtransforms.offset_copy(mtransforms.blended_transform_factory(ax.transAxes, ax.transData), fig=fig, x=-8, y=0, units="points")
species_trans = mtransforms.offset_copy(mtransforms.blended_transform_factory(ax.transAxes, ax.transData), fig=fig, x=-120, y=0, units="points")
for (spp, focal), y0 in row_y.items():
    ax.text(0, y0 + TRACK / 2, "founder" if focal == FOUNDER else focal, transform=row_trans, ha="right", va="center", fontsize=FONT,
            fontweight="bold" if focal == FOUNDER else "normal")
for spp in FOCAL_SPECIES:
    ys = [row_y[(spp, f)] for f in focal_birds]
    ax.text(0, (min(ys) + max(ys) + TRACK) / 2, species_label(spp), transform=species_trans, ha="center", va="center", rotation=90, fontsize=FONT + 3)

split = [chr_dict[c] for c in chroms].index("17")
ax.set_xticks([offset[c] + chrom_lengths[c] / 2 for c in chroms[:split]] + [(offset[chroms[split]] + total) / 2])
ax.set_xticklabels([str(chr_dict[c]) for c in chroms[:split]] + ["17-28, LGE22"], fontsize=FONT)
ax.set_xlabel("Chromosome", fontsize=FONT + 4, labelpad=8)
ax.tick_params(axis="both", length=0)
ax.set_yticks([])
ax.set_xlim(0, total)
ax.set_ylim(-0.1, max(row_y.values()) + TRACK + 0.05)
for side in ["top", "right", "left"]:
    ax.spines[side].set_visible(False)
ax.set_title(f"IBD with Daphne birds along the genome (all years): the founder and the {len(focal_birds) - 1} {species_label('conirostris')} from "
             f"Española" + (f"\n{FILTER_NOTE}" if FILTER_NOTE else ""), fontsize=FONT + 3)

fig.canvas.draw()
box = ax.get_position()
add_freq_colorbar(fig, [box.x1 + 0.01, box.y0 + box.height * 0.1, 0.01, box.height * 0.8], "Reds", vmax, label="IBD frequency")
os.makedirs(f"{OUT_DIR}/plots", exist_ok=True)
fig.savefig(f"{OUT_DIR}/plots/IBD_tracks_espanola.png", dpi=300, bbox_inches="tight", pad_inches=0.05)
plt.close(fig)
