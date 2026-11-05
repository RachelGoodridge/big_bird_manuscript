import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from common import FDR, FOCAL_SPECIES, RESULTS, TRAITS, load_bim, load_chromosomes, save_figure, set_plot_theme

set_plot_theme()
COMBINE_FROM = "17"  # chromosomes from here on are too small for their own x tick label
ALT_COLORS = {"fortis": ["#206BB5", "#0A3259"], "scandens": ["#3C9E67", "#194D30"]}  # alternating chromosome colors

# association results with each window's midpoint on a genome-wide axis (chromosomes end to end in genome order)
df = pd.read_csv(RESULTS, dtype={"chrom": str})
chr_names = load_chromosomes()
chrom_lengths = load_bim().groupby("chrom")["pos"].max().rename(index=chr_names)
chroms = list(chr_names.values())
offsets = dict(zip(chroms, np.concatenate([[0], np.cumsum(chrom_lengths[chroms].to_numpy())[:-1]])))
df["cum_pos"] = (df["window_start"] + df["window_end"]) / 2 + df["chrom"].map(offsets)
df["neg_log10_p"] = -np.log10(df["p_value"])
df["significant"] = df["q_value"] < FDR

# x ticks at each chromosome's center, with the small chromosomes sharing one label
split = chroms.index(COMBINE_FROM)
ticks = [offsets[c] + chrom_lengths[c] / 2 for c in chroms[:split]] + [(offsets[chroms[split]] + offsets[chroms[-1]] + chrom_lengths[chroms[-1]]) / 2]
tick_labels = chroms[:split] + [f"{COMBINE_FROM}-28, LGE22"]
genome_length = offsets[chroms[-1]] + chrom_lengths[chroms[-1]]


def manhattan(y_col, y_label, path, show_threshold):
    """One panel per trait (rows) and species (columns): each window's value, colored by alternating chromosome."""
    fig, axes = plt.subplots(len(TRAITS), len(FOCAL_SPECIES), figsize=(30, 16))
    for col, spp in enumerate(FOCAL_SPECIES):
        palette = {c: ALT_COLORS[spp][i % 2] for i, c in enumerate(chroms)}
        for row, trait in enumerate(TRAITS):
            ax = axes[row, col]
            sub = df[(df["species"] == spp) & (df["trait"] == trait)]
            if show_threshold:
                # genome-wide significance: the largest p-value still significant at the FDR cutoff
                sns.scatterplot(data=sub, x="cum_pos", y=y_col, hue="chrom", palette=palette, linewidth=0, s=5, alpha=0.7, legend=False, ax=ax, zorder=1)
                if sub["significant"].any():
                    ax.axhline(-np.log10(sub.loc[sub["significant"], "p_value"].max()), color="darkgray", linewidth=1, zorder=2)
                ax.set_ylim(bottom=0)
            else:
                # significant windows in red, and a line at no effect
                sns.scatterplot(data=sub[~sub["significant"]], x="cum_pos", y=y_col, hue="chrom", palette=palette, linewidth=0, s=5, alpha=0.7,
                                legend=False, ax=ax, zorder=1)
                ax.scatter(sub.loc[sub["significant"], "cum_pos"], sub.loc[sub["significant"], y_col], color="red", linewidth=0, s=8, alpha=0.7, zorder=2)
                ax.axhline(0, color="darkgray", linewidth=1, zorder=3)
            ax.set_title(f"{trait} in " + r"$\it{" + spp + "}$")
            ax.set_xlabel("Chromosome")
            ax.set_ylabel(y_label)
            ax.set_xticks(ticks)
            ax.set_xticklabels(tick_labels, fontsize=9)
            ax.set_xlim(0, genome_length)
            ax.grid(axis="x", visible=False)
            sns.despine(top=True, right=True, ax=ax)
    plt.tight_layout()
    save_figure(path, fig)


manhattan("neg_log10_p", "$-log_{10}$p-value", "manhattan_pvalues.png", show_threshold=True)
manhattan("beta", "beta (effect size)", "manhattan_beta.png", show_threshold=False)
