import importlib.util
import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from common import CONTROLS, FOUNDER, MERGE_GAP_CM, MIN_LENGTH, REPO, ORIGINAL, load_focal_birds, load_master, original_min_length, read_bim

# the figures' helpers (loaded by path since the files are all named common.py)
_spec = importlib.util.spec_from_file_location("figures_common", f"{REPO}/figures_ibis/scripts/common.py")
sys.modules["figures_common"] = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sys.modules["figures_common"])
from figures_common import GROUPS, add_bird_info, format_year_axis, italicize_species_legend, order_groups, set_plot_theme, species_label  # noqa: E402

set_plot_theme()
# with --min-length, the segments and plots of the minimum length test (see control_segments.py)
TEST = "--min-length" in sys.argv[1:]
SEGMENTS = f"{MIN_LENGTH}/trim_seg_espanola.pkl" if TEST else f"{CONTROLS}/output/trim_seg_espanola.pkl"
PLOTS = f"{MIN_LENGTH}/plots" if TEST else f"{CONTROLS}/plots"
FILTER_NOTE = f"segments > {original_min_length():.1f} cM (gaps ≤ {MERGE_GAP_CM:g} cM merged)" if TEST else ""
os.makedirs(PLOTS, exist_ok=True)
master = load_master()
focal_birds = load_focal_birds()


def save(name, fig):
    """Save (and close) a figure with the same settings as figures_ibis."""
    fig.savefig(f"{PLOTS}/{name}", dpi=450, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


# fraction of SNPs each Daphne bird of the original run shares in (filtered, trimmed) IBD with each focal bird (0 if none)
segs = pd.read_pickle(SEGMENTS)
n_snps = len(read_bim(f"{CONTROLS}/plink/dataset"))
shared = segs.groupby(["focal", "partner"])["SNP_count"].sum() / n_snps
daphne = list(pd.read_pickle(f"{ORIGINAL}/output/coef_all_founder.pkl")["ind"]) + [FOUNDER]  # every Daphne bird of the run
fractions = pd.DataFrame([(focal, bird, shared.get((focal, bird), 0.0)) for focal in focal_birds for bird in daphne if bird != focal],
                         columns=["focal", "ind", "fraction_snps_ibd"])
fractions = add_bird_info(fractions, master)
print("Daphne birds sharing IBD with each focal bird (after filtering and trimming):")
print(fractions[fractions["fraction_snps_ibd"] > 0].groupby("focal").agg(birds=("ind", "size"), max_fraction=("fraction_snps_ibd", "max"))
      .reindex(focal_birds).fillna(0).to_string())


########################### like figure 2B: each Daphne bird's IBD with each focal bird ###########################
# one panel per focal bird (the founder's is figure 2B), on the same axes as figure 2B. Each title counts the birds plotted (fortis, scandens,
# A line, B line and Big Birds of the run with known years, as figure 2B) and those sharing IBD with the focal bird after filtering and trimming
COLUMNS = 4
rows = -(-(len(focal_birds) + 1) // COLUMNS)  # one more slot for the legend
fig, axes = plt.subplots(rows, COLUMNS, figsize=(5 * COLUMNS, 4.2 * rows), sharex=True, sharey=True)
axes = axes.ravel()
for ax, focal in zip(axes, focal_birds):
    birds = order_groups(fractions[fractions["focal"] == focal], GROUPS)
    sns.scatterplot(x="Year", y="fraction_snps_ibd", hue="Species", palette=sys.modules["ibis_common"].COLOR_MAP, data=birds, s=25, alpha=0.75,
                    ax=ax, zorder=2, legend=focal == FOUNDER)
    n_ibd = (birds["fraction_snps_ibd"] > 0).sum()
    ax.set_title(f"{n_ibd} of {len(birds)} birds in IBD with {'the founder' if focal == FOUNDER else focal}", fontsize=13)
    ax.set_ylim(-0.02, 1.0)
    format_year_axis(ax)
    ax.tick_params(axis="x", rotation=30)
    ax.set_ylabel("Proportion of SNPs in IBD")
    sns.despine(top=True, right=True, ax=ax)
handles, labels = axes[0].get_legend_handles_labels()
axes[0].get_legend().remove()
for ax in axes[len(focal_birds):]:
    ax.axis("off")
legend = axes[len(focal_birds)].legend(handles, labels, loc="center", frameon=False, fontsize=13, markerscale=1.5)
italicize_species_legend(legend)
if FILTER_NOTE:
    fig.suptitle(FILTER_NOTE, fontsize=15)
fig.tight_layout()
save("year_vs_frac_ibd_espanola.png", fig)



########################### each Daphne bird's IBD with the founder vs. with each Española conirostris ###########################
# one panel per Española conirostris (2 rows of 5), with the same birds as above; the dashed line is y = x (sharing as much IBD with the
# conirostris as with the founder). The y axis is shared by all panels and zoomed to the data, as IBD with the conirostris is much lower
espanola = focal_birds[1:]
with_founder = fractions[fractions["focal"] == FOUNDER].set_index("ind")["fraction_snps_ibd"].rename("with_founder")
paired = order_groups(fractions[fractions["focal"].isin(espanola)].join(with_founder, on="ind").dropna(subset=["with_founder"]), GROUPS)
y_max = paired["fraction_snps_ibd"].max()

fig, axes = plt.subplots(2, 5, figsize=(5 * 5, 4.4 * 2))  # not shared, so every panel keeps its own y label (the limits are the same)
axes = axes.ravel()
for k, (ax, focal) in enumerate(zip(axes, espanola)):
    birds = paired[paired["focal"] == focal]
    ax.plot([0, 1], [0, 1], color="#999999", linestyle="--", linewidth=1, zorder=1)
    sns.scatterplot(x="with_founder", y="fraction_snps_ibd", hue="Species", palette=sys.modules["ibis_common"].COLOR_MAP, data=birds, s=25,
                    alpha=0.75, ax=ax, zorder=2, legend=focal == espanola[0])
    ax.set_xlim(-0.02, 1.0)
    ax.set_ylim(-y_max * 0.03, y_max * 1.08)
    ax.set_xlabel("Proportion in IBD with founder" if k >= 5 else None)  # bottom row only
    ax.set_ylabel(f"Proportion in IBD with {focal}")
    sns.despine(top=True, right=True, ax=ax)
# the legend in the empty upper part of the first panel
handles, labels = axes[0].get_legend_handles_labels()
legend = axes[0].legend(handles, labels, loc="upper left", bbox_to_anchor=(0.12, 1.0), frameon=False, fontsize=12, markerscale=1.3)
italicize_species_legend(legend)
if FILTER_NOTE:
    fig.suptitle(FILTER_NOTE, fontsize=15)
fig.tight_layout()
save("frac_ibd_founder_vs_espanola.png", fig)
