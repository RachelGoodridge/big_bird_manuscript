import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from common import (COLOR_MAP, FOCAL_SPECIES, GROUPS, IBD_LABEL, add_bill_pcs,
                    founder_ibd_per_bird, load_master, order_groups, save_figure, set_plot_theme, style_legend)

set_plot_theme()

# birds carrying IBD with the founder, with beak measurements
master = load_master()
birds = founder_ibd_per_bird(master)
birds = birds[birds["fraction_snps_ibd"] > 0].sort_values(by="ind")
birds = order_groups(birds, GROUPS)
for col in ["bill.length", "bill.depth", "bill.width"]:
    birds[col] = birds["ind"].map(master.set_index("meaningful.unique")[col])
birds = add_bill_pcs(birds.dropna(subset=["bill.length", "bill.depth", "bill.width"]))


def format_pvalue(p):
    if p < 0.001:
        return "p < 0.001"
    elif p < 0.01:
        return "p < 0.01"
    else:
        return f"p = {p:.2f}"


def plot_panel(ax, y_col, show_legend):
    sns.scatterplot(x="fraction_snps_ibd", y=y_col, hue="Species", palette=COLOR_MAP, data=birds, s=25, alpha=0.75, ax=ax, zorder=2, legend=show_legend)

    # linear fit and stats for fortis and scandens
    for y_text, (species, group) in zip([1.11, 1.03], birds[birds["Species"].isin(FOCAL_SPECIES)].groupby("Species", observed=True)):
        fit = stats.linregress(group["fraction_snps_ibd"], group[y_col])
        x_sorted = np.sort(group["fraction_snps_ibd"])
        ax.plot(x_sorted, fit.intercept + fit.slope * x_sorted, color=COLOR_MAP[species], linewidth=2)
        label = rf"$\it{{{species}}}$: $R^2$ = {fit.rvalue ** 2:.3f}, {format_pvalue(fit.pvalue)}"
        ax.text(0.03, y_text, label, transform=ax.transAxes, va="bottom", ha="left", fontsize="small", color=COLOR_MAP[species])

    ax.set_xlabel(IBD_LABEL)
    ax.set_ylabel(y_col)
    ax.set_xlim(-0.02, 0.98)
    sns.despine(top=True, right=True, left=True, bottom=True, ax=ax)


# plot fraction of SNPs against beak size and beak shape
fig, axes = plt.subplots(1, 2, figsize=(10.6, 5))
plot_panel(axes[0], "PC(Beak Size)", show_legend=False)
plot_panel(axes[1], "PC(Beak Pointedness)", show_legend=True)
style_legend(axes[1], loc="upper right")

plt.tight_layout(rect=[0, 0, 1, 0.88])
save_figure("main/beak_shape_vs_frac_ibd.png")
