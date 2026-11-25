import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from illustration_common import COLOR_MAP, save_svg, set_theme, tint

# simplified pedigree of the Big Bird lineage (on its own; figure 1B is the timeline): the founder and fortis produced the A line and the B
# line, and the B line became the Big Birds, whose fate is unknown. All positions are in the same arbitrary units and follow from the few settings below, so everything stays aligned

# ---------------- settings you may want to change ----------------
RADIUS = 1.5                 # of the circles (large enough to hold the labels)
PARENT_X = 3.3               # fortis and the founder sit this far left and right of the middle
OFFSPRING_X = 2.0            # the A line and the B line sit this far left and right of the middle
SEGMENT = 1.3                # length of every visible vertical line: between two circles, from the sibship bar to a circle, and the dashed line
DROP = RADIUS + SEGMENT      # from the line joining the parents down to the sibship bar (the bar is SEGMENT below the parents' circles)
DASH = (1.4, 1.0)            # dash and gap lengths of that arrow (in line widths)
QUESTION_GAP = 0.8           # from the arrow's tip to the middle of the question mark
BIG_BIRDS_TINT = 0.5         # the Big Birds' purple mixed with this much white, so that its black label is readable
# ------------------------------------------------------------------

# each circle: center (x, y), color, and whether the label is a species name (italic)
OFFSPRING_Y = -(DROP + SEGMENT + RADIUS)
BIG_BIRDS_Y = OFFSPRING_Y - 2 * RADIUS - SEGMENT
CIRCLES = {
    "G. fortis": dict(center=(-PARENT_X, 0.0),            color=COLOR_MAP["fortis"],  italic=True),
    "founder":   dict(center=(PARENT_X, 0.0),             color=COLOR_MAP["founder"], italic=False),
    "A line":    dict(center=(-OFFSPRING_X, OFFSPRING_Y), color=COLOR_MAP["A line"],  italic=False),
    "B line":    dict(center=(OFFSPRING_X, OFFSPRING_Y),  color=COLOR_MAP["B line"],  italic=False),
    "Big\nBirds": dict(center=(OFFSPRING_X, BIG_BIRDS_Y), color=tint(COLOR_MAP["Big Birds"], BIG_BIRDS_TINT), italic=False),
}


def draw_pedigree(ax, label_size=15, line_width=3.5):
    """The whole pedigree on one axis (equal scale on both axes, no frame)."""
    line = dict(color="black", linewidth=line_width, solid_capstyle="butt", zorder=1)
    ax.plot([-PARENT_X, PARENT_X], [0, 0], **line)                                                        # the pair
    ax.plot([0, 0], [0, -DROP], **line)                                                                   # down to their offspring
    ax.plot([-OFFSPRING_X, OFFSPRING_X], [-DROP, -DROP], **{**line, "solid_capstyle": "projecting"})      # sibship bar
    ax.plot([-OFFSPRING_X, -OFFSPRING_X], [-DROP, OFFSPRING_Y], **line)
    ax.plot([OFFSPRING_X, OFFSPRING_X], [-DROP, BIG_BIRDS_Y], **line)                                     # through the B line to the Big Birds

    # dashed arrow from the Big Birds to a question mark
    top = BIG_BIRDS_Y - RADIUS
    head_length, head_half = 0.6, 0.28
    tip = top - SEGMENT - head_length  # the dashed line is as long as the other lines, the arrow's head comes after it
    ax.plot([OFFSPRING_X, OFFSPRING_X], [top, tip + head_length], **{**line, "linestyle": (0, DASH)})
    ax.add_patch(Polygon([(OFFSPRING_X, tip), (OFFSPRING_X - head_half, tip + head_length), (OFFSPRING_X + head_half, tip + head_length)],
                         closed=True, facecolor="black", edgecolor="none", zorder=1))
    ax.text(OFFSPRING_X, tip - QUESTION_GAP, "?", ha="center", va="center_baseline", fontsize=label_size * 2.2, fontweight="bold", color="black")

    for name, c in CIRCLES.items():
        ax.add_patch(Circle(c["center"], RADIUS, facecolor=c["color"], edgecolor="none", zorder=2))
        ax.text(*c["center"], name, linespacing=1.05, multialignment="center", ha="center", va="center", fontsize=label_size, style="italic" if c["italic"] else "normal",
                color="black", zorder=3)

    ax.set_xlim(-PARENT_X - RADIUS - 0.1, PARENT_X + RADIUS + 0.1)
    ax.set_ylim(tip - QUESTION_GAP - 0.9, RADIUS + 0.1)
    ax.set_aspect("equal")
    ax.axis("off")


if __name__ == "__main__":
    set_theme()
    fig, ax = plt.subplots(figsize=(4.2, 7.5))
    draw_pedigree(ax)
    save_svg("founder_pedigree.svg", fig)
