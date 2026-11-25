import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from illustration_common import save_svg, set_theme
from admixture_panel import draw_admixture
from make_galapagos_map import EXTENT, FIG_WIDTH, draw_map
from timeline import SIZE as TIMELINE_SIZE, draw_timeline

# figure 1: (A) map of the founder's dispersal to Daphne Major, (B) timeline of the Big Bird lineage (timeline.py), (C) ADMIXTURE ancestry of
# each population

# ---------------- settings you may want to change ----------------
WIDTH = 16                 # figure width (inches)
ADMIXTURE_HEIGHT = 1.5     # height of the ADMIXTURE bars (inches)
LABEL_SIZE = 15            # text in the timeline and ADMIXTURE panels (the map's text scales with its width)
PANEL_LETTERS = True       # bold A / B / C at each panel's top left
LETTER_GAP = 0.12          # between a letter and its panel's left edge (inches; B is right at the timeline's edge, as its years start inside it)
COLUMN_GAP = 0.6           # between the map and the timeline (inches)
ROW_GAP = 1.15             # between the map and the ADMIXTURE bars (inches; holds the map's longitude labels and the ADMIXTURE legend)
# ------------------------------------------------------------------

set_theme()
# the timeline's column is exactly the width the timeline is drawn for (so its text fits) and the map takes the rest: from the figure's default
# side margins (the panels span 0.775 of the width) and the gap, as GridSpec divides them
panel_widths = WIDTH * 0.775 / (2 + COLUMN_GAP / (WIDTH / 2)) * 2
MAP_SHARE = 1 - TIMELINE_SIZE[0] / panel_widths
x0, x1, y0, y1 = EXTENT
map_width = WIDTH * MAP_SHARE
map_height = map_width * (y1 - y0) / (x1 - x0)
fig = plt.figure(figsize=(WIDTH, map_height + ADMIXTURE_HEIGHT + ROW_GAP))
# the gaps are given to GridSpec as fractions of the average panel width / height
gs = GridSpec(2, 2, figure=fig, width_ratios=[MAP_SHARE, 1 - MAP_SHARE], height_ratios=[map_height, ADMIXTURE_HEIGHT],
              wspace=COLUMN_GAP / (WIDTH / 2), hspace=ROW_GAP / ((map_height + ADMIXTURE_HEIGHT) / 2))
ax_map, ax_timeline, ax_admixture = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

draw_map(ax_map, scale=map_width / FIG_WIDTH)
# the timeline as tall as the map (known once the map is laid out), in its column
fig.canvas.draw()
box_map, box_column = ax_map.get_position(), ax_timeline.get_position()
ax_timeline.set_position([box_column.x0, box_map.y0, box_column.width, box_map.height])
timeline_size = (box_column.width * fig.get_figwidth(), box_map.height * fig.get_figheight())
print(f"timeline panel: {timeline_size[0]:.2f} x {timeline_size[1]:.2f} in (drawn for a width of {TIMELINE_SIZE[0]} in)")
draw_timeline(ax_timeline, label_size=LABEL_SIZE, size=timeline_size)
draw_admixture(ax_admixture, label_size=LABEL_SIZE)

if PANEL_LETTERS:
    # each letter just left of its panel's top left corner (A and B level with the map's top, C with the ADMIXTURE legend); positions in
    # fractions of the figure
    fig.canvas.draw()  # so the panels' final positions are known
    box_map, box_timeline, box_admixture = (ax.get_position() for ax in [ax_map, ax_timeline, ax_admixture])
    gap = LETTER_GAP / WIDTH
    legend_y = box_admixture.y1 + 0.33 / fig.get_figheight()
    letters = [("A", box_map.x0 - gap, box_map.y1, "top"), ("B", box_timeline.x0, box_map.y1, "top"),
               ("C", box_admixture.x0 - gap, legend_y, "center")]
    for letter, x, y, va in letters:
        fig.text(x, y, letter, ha="right", va=va, fontsize=24, fontweight="bold", color="black")
save_svg("figure1.svg", fig, pdf=True)
