import importlib.util
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

# shared by the hand-made figure panels (map, pedigree, deme diagram) and the figures assembled from them. Paths are found from this file's
# location (illustrations/scripts/), so scripts work from any working directory; figures are saved in illustrations/figures/ and their
# inputs are in illustrations/data/. (Not named common.py like the other folders' shared modules, so that
# figures_ibis/scripts/fig_3_combined.py can import the deme diagram next to its own common.py)
ILLUSTRATIONS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(ILLUSTRATIONS)
FIGURES = f"{ILLUSTRATIONS}/figures"
DATA = f"{ILLUSTRATIONS}/data"

# the species / lineage colors and the plot theme (font) of the other figures (loaded by path since the files are all named common.py)
if "ibis_common" not in sys.modules:
    _spec = importlib.util.spec_from_file_location("ibis_common", f"{REPO}/tune_trim_ibis/scripts/common.py")
    sys.modules["ibis_common"] = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(sys.modules["ibis_common"])
from ibis_common import COLOR_MAP, set_plot_theme


def set_theme():
    """The shared plot theme, with text saved as outlines in the svg files, so that they look the same in every viewer: saved as text ("none"
    instead of "path") it stays editable, but a viewer without the DejaVu Sans font substitutes another, whose different letter widths shift
    multi-line labels and legend entries."""
    set_plot_theme()
    plt.rcParams.update({"svg.fonttype": "path"})


def tint(color, amount):
    """A color mixed with white (amount 0: the color, 1: white), as an opaque color."""
    r, g, b = to_rgb(color)
    return (r + (1 - r) * amount, g + (1 - g) * amount, b + (1 - b) * amount)


def save_svg(name, fig, pdf=False):
    """Save (and close) a figure as an svg in illustrations/figures/. pdf: also save it as a vector pdf of the same name (for the journal), with its
    text kept as text."""
    fig.savefig(f"{FIGURES}/{name}", dpi=600, bbox_inches="tight", pad_inches=0.05)  # dpi: resolution of embedded images (the photo)
    if pdf:
        with plt.rc_context({"pdf.fonttype": 42}):
            fig.savefig(f"{FIGURES}/{os.path.splitext(name)[0]}.pdf", dpi=600, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
