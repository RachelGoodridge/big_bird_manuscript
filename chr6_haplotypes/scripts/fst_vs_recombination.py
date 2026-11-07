import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec
from scipy.optimize import least_squares
from scipy.stats import spearmanr
from common import CHROM, END, FST_WINDOWS, GMAP, OUTPUT, PEAK_REGIONS, START, WORK, save_figure, set_plot_theme

RESID_Z_THRESH = 3.5  # FST peaks: windows with FST at least this many SDs above that expected from their recombination rate
LOSS_SCALE = 0.05   # residual size (FST units) beyond which the robust fit down-weights windows

set_plot_theme()

########################### LOAD DATA ###########################
windows = pd.read_csv(FST_WINDOWS, sep="\t")
gmap = pd.read_csv(GMAP, sep="\t")

# recombination rate in each FST window (cM/Mb) from the interpolated genetic map
cm_start, cm_end = (np.interp(windows[c], gmap["pos"], gmap["cM"]) for c in ["start", "end"])
windows["recomb_rate"] = (cm_end - cm_start) / ((windows["end"] - windows["start"]) / 1e6)
rho, p = spearmanr(windows["recomb_rate"], windows["fst"])
print(f"Spearman rho (recombination rate vs window FST): {rho:.3f} (p = {p:.2g})")

########################### FST EXPECTED FROM RECOMBINATION ###########################
# linked selection predicts FST rising as recombination falls: FST = a + b / (r + c), fit with a robust (soft L1) loss so the FST peaks
# themselves pull less on the expected curve; residual = observed - expected FST
def inverse_model(params, r):
    a, b, c = params
    return a + b / (r + c)

rate = windows["recomb_rate"].to_numpy()
fit = least_squares(lambda prm: inverse_model(prm, rate) - windows["fst"].to_numpy(), x0=[0.1, 0.05, 0.1], bounds=([-1, 0, 1e-3], [1, 10, 10]),
                    loss="soft_l1", f_scale=LOSS_SCALE)
a, b, c = fit.x
fit_label = rf"$F_{{ST}} = {a:.3f} + \frac{{{b:.3f}}}{{r + {c:.3f}}}$"
print(f"Fit: FST = {a:.3f} + {b:.3f} / (r + {c:.3f})")
windows["fst_expected"] = inverse_model(fit.x, rate)
windows["resid"] = windows["fst"] - windows["fst_expected"]
windows["resid_z"] = (windows["resid"] - windows["resid"].mean()) / windows["resid"].std()
windows.to_csv(f"{OUTPUT}/fst_recombination_windows.tsv", sep="\t", index=False, float_format="%.5f")


def merge_windows(df):
    """Overlapping windows merged into regions (start, end)."""
    regions = []
    for _, w in df.iterrows():
        if regions and w["start"] <= regions[-1][1]:
            regions[-1][1] = max(regions[-1][1], w["end"])
        else:
            regions.append([w["start"], w["end"]])
    return regions


# FST peaks: overlapping windows with residual Z >= RESID_Z_THRESH merged into regions (used for the genotype heatmap)
peaks = []
for s, e in merge_windows(windows[windows["resid_z"] >= RESID_Z_THRESH]):
    sub = windows[(windows["start"] >= s) & (windows["end"] <= e)]
    peaks.append({"start": s, "end": e, "max_fst": sub["fst"].max(), "max_resid_z": sub["resid_z"].max(), "mean_recomb_rate": sub["recomb_rate"].mean()})
peaks = pd.DataFrame(peaks, columns=["start", "end", "max_fst", "max_resid_z", "mean_recomb_rate"]).astype({"start": int, "end": int})
peaks.to_csv(PEAK_REGIONS, sep="\t", index=False, float_format="%.3f")
print(f"FST peaks (residual Z >= {RESID_Z_THRESH}):")
print(peaks.to_string(index=False))

########################### PLOT ###########################
# left: FST against recombination rate; right: windowed FST and FST above that expected from recombination, along the region
peak = windows["resid_z"] >= RESID_Z_THRESH
mid_mb = windows["mid"] / 1e6

fig = plt.figure(figsize=(22, 7.5))
gs = GridSpec(2, 2, figure=fig, width_ratios=[1, 2.4], wspace=0.15, hspace=0.35)
ax_scatter = fig.add_subplot(gs[:, 0])
axes = [fig.add_subplot(gs[i, 1]) for i in range(2)]

ax_scatter.scatter(windows.loc[~peak, "recomb_rate"], windows.loc[~peak, "fst"], s=14, color="#1a5276", alpha=0.5, linewidth=0, label="100 kb windows")
ax_scatter.scatter(windows.loc[peak, "recomb_rate"], windows.loc[peak, "fst"], s=26, color="#c0392b", linewidth=0, label=f"FST peaks (residual Z ≥ {RESID_Z_THRESH})")
r_grid = np.linspace(rate.min(), rate.max(), 300)
ax_scatter.plot(r_grid, inverse_model(fit.x, r_grid), color="black", linewidth=1.2, linestyle="--", label=f"Expected {fit_label}")
ax_scatter.set_xlabel("Recombination rate (cM/Mb)")
ax_scatter.set_ylabel("Hudson FST (100 kb windows)")
ax_scatter.set_xlim(left=0)
ax_scatter.set_ylim(bottom=0)
ax_scatter.set_box_aspect(1)
ax_scatter.set_title(f"Spearman ρ = {rho:.2f}", loc="left")
ax_scatter.legend(loc="upper right", fontsize=11, frameon=False)

axes[0].plot(mid_mb, windows["fst"], color="#1a5276", linewidth=1, label="100 kb windows (step 50 kb)", zorder=2)
axes[0].plot(mid_mb, windows["fst_expected"], color="black", linewidth=1, linestyle="--", label="Expected from recombination", zorder=3)
axes[0].set_ylabel("Hudson FST")
axes[0].set_title("FST, conirostris + magnirostris vs scandens (no IBD)", loc="left")
axes[0].set_ylim(0, 1)
axes[0].legend(loc="upper right", fontsize=10, frameon=False, ncol=2)

axes[1].plot(mid_mb, windows["resid_z"], color="#1a5276", linewidth=1)
axes[1].axhline(RESID_Z_THRESH, color="red", linestyle="--", linewidth=0.8)
axes[1].axhline(0, color="lightgray", linewidth=0.8, zorder=0)
axes[1].set_ylabel("Z score")
axes[1].set_title(f"FST above that expected from recombination (peaks: Z ≥ {RESID_Z_THRESH}, shaded)", loc="left")
axes[1].set_xlabel(f"Chromosome {CHROM} position (Mb)")

for ax in axes:
    for s, e in zip(peaks["start"], peaks["end"]):
        ax.axvspan(s / 1e6, e / 1e6, color="#c0392b", alpha=0.2, linewidth=0, zorder=0)
    ax.set_xlim(START / 1e6, END / 1e6)
axes[0].tick_params(labelbottom=False)
for ax in [ax_scatter] + axes:
    ax.grid(False)
    sns.despine(ax=ax)
save_figure("fst_recombination.png", fig)
