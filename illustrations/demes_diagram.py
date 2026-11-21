import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from illustration_common import COLOR_MAP, save_svg, set_theme, tint

# schematic of the populations on Daphne Major through time (figure 3): each population is a tube whose width shows its relative size, thick arrows
# are the founding and outcrossing events, and thin arrows are ongoing gene flow. Years run down the page; the horizontal positions and the widths
# are schematic (not to scale)

# ---------------- settings you may want to change ----------------
START, END = 1981, 2022
FOUNDER_YEARS = (1981, 1994)  # years the founder was on Daphne
A_LINE_START = 1987           # founder x a fortis-scandens hybrid (drawn as fortis, to keep the diagram simple)
B_LINE_START = 1991           # founder x fortis
BIG_BIRDS_START = 2006        # the B line from here on
OUTCROSS = 2011               # Big Birds x scandens
FADE_START = 2012             # the A line and the Big Birds fade from here to END
FADE_SPREAD = 0.45            # how the fade runs (a bell-curve-like rise): smaller turns white sooner after FADE_START; it starts and ends
                              # gradually, with no visible edge at FADE_START
GENE_FLOW_START = 2000        # fortis <-> scandens gene flow, more and more frequent from here on

# each tube: color, horizontal center, and half-width at the given years (straight lines in between)
TUBES = {
    "scandens":  dict(color=COLOR_MAP["scandens"],  center=1.50,  years=[START, END],                         half_width=[0.55, 0.90]),
    "fortis":    dict(color=COLOR_MAP["fortis"],    center=4.60,  years=[START, END],                         half_width=[0.95, 0.65]),
    "A line":    dict(color=COLOR_MAP["A line"],    center=7.40,  years=[A_LINE_START, 1992, END],            half_width=[0.18, 0.45, 0.15]),
    "B line":    dict(color=COLOR_MAP["B line"],    center=10.20, years=[B_LINE_START, BIG_BIRDS_START],      half_width=[0.18, 0.36]),
    "Big Birds": dict(color=COLOR_MAP["Big Birds"], center=10.20, years=[BIG_BIRDS_START, FADE_START, END],   half_width=[0.36, 0.48, 0.48]),
    "founder":   dict(color=COLOR_MAP["founder"],   center=12.90, years=list(FOUNDER_YEARS),                  half_width=[0.07, 0.07]),
}
FADING = ["A line", "Big Birds"]
FILL_TINT = 0.55     # how much lighter the tubes' fill is than their outline
THICK, THIN = 3.2, 0.7  # line widths of the events and of the gene flow
HALO_OPACITY = 0.6      # of the white band under an event's arrow where it crosses a tube (0: no band)
A_GENE_FLOW_GAP = 1.5   # years between the gene flow arrows of the A line (and of the B line)
N_SPECIES_ARROWS = 15   # fortis <-> scandens arrows, the last gaps LAST_GAP times the first
LAST_GAP = 0.2
# ------------------------------------------------------------------


def edges(name, year):
    """Left and right edge of a tube in a year."""
    t = TUBES[name]
    half = np.interp(year, t["years"], t["half_width"])
    return t["center"] - half, t["center"] + half


def fade(year):
    """How far a fading tube (and its arrows) has turned white in a year: 0 until FADE_START, 1 from END."""
    t = np.clip((np.asarray(year) - FADE_START) / (END - FADE_START), 0, 1)
    return (1 - np.exp(-(t / FADE_SPREAD) ** 2)) / (1 - np.exp(-(1 / FADE_SPREAD) ** 2))


def draw_tube(ax, name, zorder):
    """A tube: lighter fill with its two sides outlined. Fading tubes are covered from FADE_START to END by white of rising opacity (on the
    white page this looks the same as the tube turning transparent, without a seam where the fade starts)."""
    t = TUBES[name]
    years = np.array(t["years"])
    left, right = np.array([edges(name, y) for y in years]).T
    ax.add_patch(Polygon(np.column_stack([np.concatenate([left, right[::-1]]), np.concatenate([years, years[::-1]])]), closed=True,
                         facecolor=tint(t["color"], FILL_TINT), edgecolor="none", zorder=zorder))
    for side in (left, right):
        ax.plot(side, years, color=t["color"], linewidth=1.8, solid_capstyle="butt", zorder=zorder + 0.1)
    if name in FADING:
        cover = np.ones((256, 1, 4))
        cover[:, 0, 3] = fade(np.linspace(FADE_START, END + 0.2, 256))
        ax.imshow(cover, extent=[left.min() - 0.1, right.max() + 0.1, END + 0.2, FADE_START], aspect="auto", zorder=zorder + 0.2,
                  interpolation="bilinear")


def arrow(ax, x_from, x_to, year, color, width, zorder, dashed=False, halo=False, gap=0.06):
    """Horizontal arrow in a year, stopping just short of where it points (halo: a translucent white band under it, so it reads clearly where it
    crosses a tube)."""
    direction = np.sign(x_to - x_from)
    head_length, head_half = (0.30, 0.38) if width == THICK else (0.17, 0.2)
    tip = x_to - direction * gap
    shaft = [x_from + direction * gap, tip - direction * head_length * 0.9]
    if halo:
        ax.plot(shaft, [year, year], color="white", alpha=HALO_OPACITY, linewidth=width + 1.6, solid_capstyle="butt", zorder=zorder)
    ax.plot(shaft, [year, year], color=color, linewidth=width, solid_capstyle="butt", linestyle=(0, (3, 1.6)) if dashed else "-", zorder=zorder)
    ax.add_patch(Polygon([(tip, year), (tip - direction * head_length, year - head_half), (tip - direction * head_length, year + head_half)],
                         closed=True, facecolor=color, edgecolor="none", zorder=zorder))


def gene_flow(ax, deme_a, deme_b, years, zorder, fading=False):
    """Thin arrows between two neighboring tubes (deme_a on the left), alternating in direction, each the color of the deme it comes from.
    fading: arrows after FADE_START turn white as the fading tubes do (at the same rate as their white cover)."""
    for k, year in enumerate(years):
        x_a, x_b = edges(deme_a, year)[1], edges(deme_b, year)[0]
        whiteness = 0.25
        if fading:
            whiteness += (1 - whiteness) * fade(year)
        if k % 2:
            arrow(ax, x_a, x_b, year, tint(TUBES[deme_a]["color"], whiteness), THIN, zorder)
        else:
            arrow(ax, x_b, x_a, year, tint(TUBES[deme_b]["color"], whiteness), THIN, zorder)


def avoid(years, events, margin=0.45):
    """Years of gene flow arrows, without those that would sit on top of an event."""
    return [y for y in years if all(abs(y - e) > margin for e in events)]


def draw_demes(ax, label_size=15):
    """The whole diagram on one axis."""
    events = [A_LINE_START, B_LINE_START, OUTCROSS]

    for name in TUBES:
        draw_tube(ax, name, zorder=2)
    # gene flow between fortis and the A line (until END, fading with the A line), and between fortis and the B line (halfway between the A
    # line's arrows, crossing over the A line)
    gene_flow(ax, "fortis", "A line", avoid(np.arange(A_LINE_START + 1, END, A_GENE_FLOW_GAP), events), zorder=3, fading=True)
    b_years = np.arange(A_LINE_START + 1 + A_GENE_FLOW_GAP / 2, BIG_BIRDS_START, A_GENE_FLOW_GAP)
    gene_flow(ax, "fortis", "B line", [y for y in avoid(b_years, events) if y > B_LINE_START], zorder=3)

    # fortis <-> scandens: the gaps between arrows shrink at a constant rate from GENE_FLOW_START to END
    gaps = np.linspace(1, LAST_GAP, N_SPECIES_ARROWS - 1)
    years = GENE_FLOW_START + np.concatenate([[0], np.cumsum(gaps)]) * (END - 0.6 - GENE_FLOW_START) / gaps.sum()
    gene_flow(ax, "scandens", "fortis", avoid(years, events, margin=0.3), zorder=3)

    # events: fortis and the founder start the A line and then the B line; the Big Birds outcross to scandens
    for start, line in [(A_LINE_START, "A line"), (B_LINE_START, "B line")]:
        left, right = edges(line, start)
        arrow(ax, edges("fortis", start)[1], left, start, COLOR_MAP["fortis"], THICK, zorder=4, halo=line == "B line")
        arrow(ax, edges("founder", start)[0], right, start, COLOR_MAP["founder"], THICK, zorder=4)
    arrow(ax, edges("Big Birds", OUTCROSS)[0], edges("scandens", OUTCROSS)[1], OUTCROSS, COLOR_MAP["Big Birds"], THICK, zorder=4, dashed=True,
          halo=True)

    # names above each tube (the Big Birds beside theirs)
    for name in ["scandens", "fortis", "A line", "B line", "founder"]:
        t = TUBES[name]
        species = name in ["scandens", "fortis"]  # species names in italics, with the genus
        ax.text(t["center"], t["years"][0] - 0.5, f"G. {name}" if species else name, ha="center", va="bottom", fontsize=label_size,
                style="italic" if species else "normal")
    ax.text(edges("Big Birds", 2008.5)[1] + 0.2, 2008.5, "Big Birds", ha="left", va="center", fontsize=label_size)

    # year axis on the left
    ax.set_ylim(END + 0.3, START - 2)
    ax.set_xlim(0.2, 13.2)  # just past the founder (its name may stick out; the saved figure includes it)
    ax.set_yticks([START, A_LINE_START, B_LINE_START, BIG_BIRDS_START, OUTCROSS, END])
    ax.set_xticks([])
    ax.grid(False)
    for side, spine in ax.spines.items():
        spine.set_visible(side == "left")
    ax.spines["left"].set_bounds(START, END)
    ax.spines["left"].set_color("black")
    ax.tick_params(axis="y", left=True, length=6, width=1.2, color="black")


if __name__ == "__main__":
    set_theme()
    fig, ax = plt.subplots(figsize=(10, 8))
    draw_demes(ax)
    save_svg("figure3_demes_diagram.svg", fig)
