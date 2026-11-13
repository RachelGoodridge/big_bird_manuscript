import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from illustration_common import save_svg, set_theme
from admixture_panel import draw_admixture
from make_galapagos_map import EXTENT, FIG_WIDTH, draw_map
from pedigree import draw_pedigree

# figure 1: (A) map of the founder's dispersal to Daphne Major, (B) pedigree of the Big Bird lineage, (C) ADMIXTURE ancestry of each population

# ---------------- settings you may want to change ----------------
WIDTH = 16                 # figure width (inches)
MAP_SHARE = 0.70           # share of the width taken by the map (the pedigree takes the rest)
ADMIXTURE_HEIGHT = 1.5     # height of the ADMIXTURE bars (inches)
LABEL_SIZE = 15            # text in the pedigree and ADMIXTURE panels (the map's text scales with its width)
PANEL_LETTERS = True       # bold A / B / C at each panel's top left
LETTER_GAP = 0.12          # between a letter and its panel's left edge (inches)
COLUMN_GAP = 0.6           # between the map and the pedigree (inches)
ROW_GAP = 1.15             # between the map and the ADMIXTURE bars (inches; holds the map's longitude labels and the ADMIXTURE legend)
# ------------------------------------------------------------------

set_theme()
x0, x1, y0, y1 = EXTENT
map_width = WIDTH * MAP_SHARE
map_height = map_width * (y1 - y0) / (x1 - x0)
fig = plt.figure(figsize=(WIDTH, map_height + ADMIXTURE_HEIGHT + ROW_GAP))
# the gaps are given to GridSpec as fractions of the average panel width / height
gs = GridSpec(2, 2, figure=fig, width_ratios=[MAP_SHARE, 1 - MAP_SHARE], height_ratios=[map_height, ADMIXTURE_HEIGHT],
              wspace=COLUMN_GAP / (WIDTH / 2), hspace=ROW_GAP / ((map_height + ADMIXTURE_HEIGHT) / 2))
ax_map, ax_pedigree, ax_admixture = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])

draw_map(ax_map, scale=map_width / FIG_WIDTH)
draw_pedigree(ax_pedigree, label_size=LABEL_SIZE)
ax_pedigree.set_anchor("N")  # the pedigree's top level with the map's
draw_admixture(ax_admixture, label_size=LABEL_SIZE)

if PANEL_LETTERS:
    # each letter just left of its panel's top left corner (A and B level with the map's top, C with the ADMIXTURE legend); positions in
    # fractions of the figure
    fig.canvas.draw()  # so the panels' final positions are known
    box_map, box_pedigree, box_admixture = (ax.get_position() for ax in [ax_map, ax_pedigree, ax_admixture])
    gap = LETTER_GAP / WIDTH
    legend_y = box_admixture.y1 + 0.33 / fig.get_figheight()
    letters = [("A", box_map.x0 - gap, box_map.y1, "top"), ("B", box_pedigree.x0 - gap, box_map.y1, "top"),
               ("C", box_admixture.x0 - gap, legend_y, "center")]
    for letter, x, y, va in letters:
        fig.text(x, y, letter, ha="right", va=va, fontsize=24, fontweight="bold", color="black")
save_svg("figure1.svg", fig, pdf=True)
