import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyBboxPatch, Polygon
from matplotlib.textpath import TextPath
from illustration_common import COLOR_MAP, save_svg, set_theme, tint
from ibis_common import species_label
from demes_diagram import A_LINE_START, B_LINE_START, BIG_BIRDS_START, END, FADE_START, FOUNDER_YEARS

# timeline of the Big Bird lineage (an alternative figure 1B, a simplified figure 3C): the founder started the A line (with a hybrid; later absorbed
# into G. fortis) and the B line (with a G. fortis), and the B line became the reproductively isolated Big Birds, whose fate is unknown.
# Time runs downward. Horizontal positions are in inches, vertical ones in years

# ---------------- settings you may want to change ----------------
SIZE = (3.7, 6.0)             # figure width, height (inches); the width just holds the three bars at the text sizes below
YEARS = (1978.8, 2024.8)      # top and bottom of the drawing
# each lineage's first and last year, all taken from the settings of figure 3C (demes_diagram.py) so that the two figures agree: the A line and
# the Big Birds end where they start to fade there, and the question mark is at that diagram's last year, when no Big Birds were seen
FOUNDER = FOUNDER_YEARS
A_LINE = (A_LINE_START, FADE_START)
B_LINE = (B_LINE_START, BIG_BIRDS_START)
BIG_BIRDS = (BIG_BIRDS_START, FADE_START)
SURVEY = END
AXIS_X = 0.66                 # the time axis
BAR_X = {"A line": 1.27, "founder": 2.25, "B line": 3.23}   # middle of each lineage's bar
BAR_WIDTH = 0.62              # of the A line, the B line and the Big Birds
FOUNDER_WIDTH = 0.40          # of the founder (a single bird)
CORNER = 0.13                 # radius of the bars' rounded corners
GAP = 0.07                    # between a shape and the text next to it
BIG_BIRDS_TINT = 0.5          # the Big Birds' purple mixed with this much white, as in the pedigree, so that its black label is readable
# ------------------------------------------------------------------

def draw_timeline(ax, label_size=15, size=SIZE):
    """The whole timeline on one axis that is size inches (width, height): the width should stay close to SIZE's, as the horizontal positions
    are in inches, but any height works. The lineages' names are a little smaller than label_size
    (upright, inside their bars; the founder's above its bar, which is too narrow for it), the years and notes smaller again."""
    name_size, year_size, note_size = label_size * 0.85, label_size * 0.8, label_size * 0.7
    years_per_inch = (YEARS[1] - YEARS[0]) / size[1]
    pt = years_per_inch / 72  # one point, in years
    half = BAR_WIDTH / 2
    a, founder, b = BAR_X["A line"], BAR_X["founder"], BAR_X["B line"]
    big_birds_color = tint(COLOR_MAP["Big Birds"], BIG_BIRDS_TINT)
    name = dict(ha="center", fontsize=name_size, color="black", linespacing=1.05, multialignment="center", zorder=3)
    note = dict(fontsize=note_size, color="black", linespacing=1.1)

    def bar(x, first, last, color, label=None, width=BAR_WIDTH):
        ax.add_patch(FancyBboxPatch((x - width / 2, first), width, last - first, boxstyle=f"round,pad=0,rounding_size={CORNER}",
                                    mutation_aspect=years_per_inch, facecolor=color, edgecolor="none", zorder=2))
        if label:
            ax.text(x, (first + last) / 2, label, va="center", **name)

    def end_bar(x, year):
        """Black line across a bar's end."""
        ax.plot([x - half - 0.08, x + half + 0.08], [year, year], color="black", lw=2, solid_capstyle="butt", zorder=4)

    # time axis, with a tick at each event
    ax.annotate("", xy=(AXIS_X, YEARS[1] - 0.6), xytext=(AXIS_X, YEARS[0] + 0.6),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.5, shrinkA=0, shrinkB=0, mutation_scale=14))
    for year in sorted({FOUNDER[0], A_LINE[0], B_LINE[0], BIG_BIRDS[0], BIG_BIRDS[1], SURVEY}):
        ax.plot([AXIS_X - 0.06, AXIS_X + 0.06], [year, year], color="black", lw=1.5, solid_capstyle="butt")
        ax.text(AXIS_X - 0.06 - GAP, year, str(year), ha="right", va="center", fontsize=year_size, color="black")

    # lineages
    bar(founder, *FOUNDER, COLOR_MAP["founder"], width=FOUNDER_WIDTH)
    ax.text(founder, FOUNDER[0] - 4 * pt, "founder", va="bottom", **name)
    bar(a, *A_LINE, COLOR_MAP["A line"], "A line")
    bar(b, *B_LINE, COLOR_MAP["B line"], "B line")
    bar(b, *BIG_BIRDS, big_birds_color, "Big\nBirds")

    # an arrow from the founder to each lineage it started, in the lineage's first year, ending where the bar's corner would be if it were
    # not rounded (so just short of the bar)
    arrow = dict(arrowstyle="-|>", color="black", lw=1.5, shrinkA=0, shrinkB=0, mutation_scale=12)
    ax.annotate("", xy=(a + half, A_LINE[0]), xytext=(founder - FOUNDER_WIDTH / 2, A_LINE[0]), arrowprops=arrow, zorder=4)
    ax.annotate("", xy=(b - half, B_LINE[0]), xytext=(founder + FOUNDER_WIDTH / 2, B_LINE[0]), arrowprops=arrow, zorder=4)

    # black lines across the lineages' ends: the A line's (absorbed into fortis), the B line's (reproductive isolation: the Big Birds from here
    # on) and the Big Birds' (last observed)
    end_bar(a, A_LINE[1])
    ax.text(a, A_LINE[1] + 6 * pt, "absorbed into\n" + species_label("fortis"), ha="center", va="top", multialignment="center", **note)
    end_bar(b, BIG_BIRDS[0])
    ax.text(b - half - 0.08 - GAP, BIG_BIRDS[0], "reproductive\nisolation", ha="right", va="center", multialignment="right", **note)  # left of the line
    end_bar(b, BIG_BIRDS[1])
    ax.text(b - GAP, BIG_BIRDS[1] + 6 * pt, "last\nobserved", ha="right", va="top", multialignment="right", **note)  # under the line, left of the dashes

    # unknown fate: a dashed arrow from the Big Birds' end to a question mark at the survey in which none were seen. The arrow starts as far
    # below the black line as it ends above the question mark (measured from the glyph's outline, which is centered on the survey year)
    question_size, gap, line_half = label_size * 2.2, 3 * pt, 1 * pt  # line_half: half the black line's width
    glyph = TextPath((0, 0), "?", size=question_size, prop=FontProperties(family=plt.rcParams["font.family"], weight="bold")).get_extents()
    ax.text(b, SURVEY + (glyph.y0 + glyph.y1) / 2 * pt, "?", ha="center", va="baseline", fontsize=question_size, fontweight="bold", color="black")
    head_length, head_half = 0.16 * years_per_inch, 0.085
    tip = SURVEY - glyph.height / 2 * pt - gap
    ax.plot([b, b], [BIG_BIRDS[1] + line_half + gap, tip - head_length], color=big_birds_color, lw=3, linestyle=(0, (1.4, 1.0)), solid_capstyle="butt", zorder=1)
    ax.add_patch(Polygon([(b, tip), (b - head_half, tip - head_length), (b + head_half, tip - head_length)], closed=True,
                         facecolor=big_birds_color, edgecolor="none", zorder=1))

    ax.set_xlim(0, size[0])
    ax.set_ylim(YEARS[1], YEARS[0])  # time runs downward
    ax.axis("off")


if __name__ == "__main__":
    set_theme()
    fig = plt.figure(figsize=SIZE)
    draw_timeline(fig.add_axes([0, 0, 1, 1]))
    save_svg("figure1_timeline.svg", fig)
