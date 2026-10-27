import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import seaborn as sns
from common import DATA, FOCAL_SPECIES, GENOTYPING, ibd_window_frequencies, load_chromosomes, load_master, save_figure, set_plot_theme, style_legend

set_plot_theme()

# proportion of fortis / scandens alive each year in IBD with the founder in each genomic window
freqs, windows = ibd_window_frequencies(load_master())

# windows within each of the 6 large-effect beak loci
loci_pos = pd.read_csv(f"{GENOTYPING}/data/loci_positions.csv")
chrom_code = {f"chr{name}": code for code, name in load_chromosomes().items()}
loci_wins = {}
for _, locus in loci_pos.iterrows():
    inside = (windows["chrom"] == chrom_code[locus["chrom"]]) & (windows["start"] > locus["pos_start"]) & (windows["end"] < locus["pos_end"])
    loci_wins[locus["locus"].split("_")[-1]] = windows.loc[inside, "name"].to_list()
locus_windows = [w for wins in loci_wins.values() for w in wins]
colors = dict(zip(loci_wins, plt.rcParams["axes.prop_cycle"].by_key()["color"]))
print({locus: len(wins) for locus, wins in loci_wins.items()})


# one line per window over time
def window_lines(freq, cols, **kwargs):
    return LineCollection([np.column_stack([freq.index, freq[c].to_numpy()]) for c in cols], **kwargs)


# every window (gray), windows in the beak loci (colors), and the mean across all windows (black)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
for ax, spp in zip(axes, FOCAL_SPECIES):
    freq = freqs[spp]
    ax.add_collection(window_lines(freq, freq.columns[~freq.columns.isin(locus_windows)], color="lightgray", linewidth=0.5, alpha=0.25, zorder=1))
    for locus, wins in loci_wins.items():
        ax.add_collection(window_lines(freq, wins, color=colors[locus], linewidth=0.75, alpha=0.6, zorder=2))
    ax.plot(freq.index, freq.mean(axis=1), color="black", linewidth=2.0, alpha=0.9, zorder=3)

    ax.set_title(spp, fontstyle="italic")
    ax.set_xlabel(None)
    ax.set_xticks(np.arange(1990, 2025, 5))
    ax.set_xlim(1988, 2024)
    ax.set_ylim(0, max(f.to_numpy().max() for f in freqs.values()) * 1.02)
    sns.despine(top=True, right=True, ax=ax)

axes[0].set_ylabel("Proportion of Individuals in IBD")
style_legend(axes[0], handles=[Line2D([0], [0], color=colors[locus], lw=0.75, label=locus) for locus in loci_wins], loc="upper left")
plt.tight_layout()
save_figure("supplement/IBDrw_allele_freqs_over_time.png", fig)
