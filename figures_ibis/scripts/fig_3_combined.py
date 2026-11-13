import sys
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from common import REPO, add_panel_letters, load_master, load_max_kin, plot_max_kin, save_figure, set_plot_theme

# the schematic of the populations on Daphne Major through time (drawn by illustrations/demes_diagram.py)
sys.path.append(f"{REPO}/illustrations")
from demes_diagram import draw_demes

# figure 3: (A, B) kinship of every fortis / scandens to its closest relative among the founder's descendants,
# (C) schematic of the populations and the gene flow between them

# ---------------- settings you may want to change ----------------
KINSHIP_SIZE = (12, 5)      # size the two kinship panels are laid out for on their own; their proportions and spacing are kept (inches)
DIAGRAM_SIZE = (10, 8)      # size of the schematic's own figure (illustrations/demes_diagram.py), whose proportions are kept
SCALE = 1.3                 # everything in the kinship plot (panels, text, points) is this many times larger than on its own
DIAGRAM_GAP = 1.15          # blank space between the kinship panels and the schematic (inches; holds the degrees of relatedness)
DIAGRAM_LABEL_SIZE = 15     # names in the schematic
# ------------------------------------------------------------------

set_plot_theme()
pivot, cutoffs = load_max_kin(load_master())


def kinship_axes(fig):
    """The two kinship panels on a figure of their own size (KINSHIP_SIZE), with matplotlib's tight layout."""
    axes = fig.subplots(ncols=2, sharex=True, sharey=True)
    plot_max_kin(axes, pivot, cutoffs)
    fig.tight_layout()
    return axes


# where the two kinship panels sit in the plot on its own, and the plot area of the schematic's own figure (fractions of those figures)
probe = plt.figure(figsize=KINSHIP_SIZE)
kinship_boxes = [ax.get_position() for ax in kinship_axes(probe)]
plt.close(probe)
probe = plt.figure(figsize=DIAGRAM_SIZE)
diagram_box = probe.add_subplot().get_position()
plt.close(probe)

# the kinship plot, enlarged, takes the left of the figure; the schematic is scaled to the height of the kinship panels
kinship_width, height = KINSHIP_SIZE[0] * SCALE, KINSHIP_SIZE[1] * SCALE
panel_height = kinship_boxes[0].height * height
diagram_width = panel_height * (diagram_box.width * DIAGRAM_SIZE[0]) / (diagram_box.height * DIAGRAM_SIZE[1])
width = kinship_width + DIAGRAM_GAP + diagram_width
plt.rcParams.update({key: plt.rcParams[key] * SCALE for key in ["font.size", "axes.labelsize", "axes.titlesize", "xtick.labelsize", "ytick.labelsize",
                                                                "legend.fontsize", "lines.markersize"]})
fig = plt.figure(figsize=(width, height))
ax_fortis = fig.add_axes([kinship_boxes[0].x0 * kinship_width / width, kinship_boxes[0].y0, kinship_boxes[0].width * kinship_width / width,
                          kinship_boxes[0].height])
ax_scandens = fig.add_axes([kinship_boxes[1].x0 * kinship_width / width, kinship_boxes[1].y0, kinship_boxes[1].width * kinship_width / width,
                            kinship_boxes[1].height], sharex=ax_fortis, sharey=ax_fortis)
ax_demes = fig.add_axes([(kinship_width + DIAGRAM_GAP) / width, kinship_boxes[0].y0, diagram_width / width, kinship_boxes[0].height])

plot_max_kin([ax_fortis, ax_scandens], pivot, cutoffs, point_size=50 * SCALE ** 2)
ax_scandens.tick_params(labelleft=False)  # the two kinship panels share the y axis
ax_scandens.set_ylabel(None)
draw_demes(ax_demes, label_size=DIAGRAM_LABEL_SIZE)

add_panel_letters(fig, [ax_fortis, ax_scandens, ax_demes], fontsize=28, pad=14, shift={1: 0.3})  # B a little left of its y axis
save_figure("main/figure3.png", fig, pdf=True)
