"""
Galápagos map showing the dispersal of the G. conirostris founder
from Española to Daphne Major.

Setup (once):  pip install geopandas matplotlib
Run:           python make_galapagos_map.py
Needs:         galapagos_land.geojson and conirostris_head.png
Writes:        galapagos_map.svg
Also used by:  fig_1_combined.py, which draws the map as panel A with draw_map()

Coastlines: derived from OpenStreetMap land polygons
            (© OpenStreetMap contributors, ODbL) via the geo-maps project.
"""
import numpy as np
import geopandas as gpd
import os
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from pyproj import Geod
from shapely.geometry import Polygon, LineString, Point
from illustration_common import set_theme

# ---------------- settings you may want to change ----------------
SHOW_TOWNS  = False              # True adds Puerto Ayora etc.
EXTENT      = (-91.74, -89.14, -1.50, 0.21)   # lon_min, lon_max, lat_min, lat_max (0.74 in the north would add Pinta, Marchena and Genovesa)
FIG_WIDTH   = 7.0                # inches (~180 mm, a typical double column); text and line sizes below are for this width
LAND, COAST = "#a5c98b", "#6f9f63"            # greens, as in the original map
OCEAN       = "#a8d8f0"                       # light blue
TEXT        = "#2b2b2b"
SEA_TEXT    = "#4f7f9c"          # "Pacific Ocean" and the title
TITLE       = "Islas Galápagos"  # top right corner ("" for none)
RED         = "#e8212b"          # route colour
DIST_LABEL  = None               # text along the route; None: the length of the route as drawn (the red dashed line), in km
COAST_GAP_KM = 0.5               # how far offshore the route runs along Santa Cruz
ROUTE_LW    = 1.2                # route line width (points)
ROUTE_DASH  = (2.2, 1.6)         # dash length, gap length (in line widths)
ARROW_HEAD  = (0.42, 0.2)        # arrowhead length, half-width
ESPANOLA    = (-89.677, -1.377)  # lon, lat
DAPHNE      = (-90.372, -0.422)
ROUTE_ENTRY = (-90.20, -0.735)   # where the route reaches the Santa Cruz coast (SE corner)
ROUTE_EXIT  = (-90.325, -0.455)    # where it leaves the coast for Daphne Major
BIRD_IMAGE  = "conirostris_head.png"
BIRD_BOX    = (0.845, 0.095, 0.13)   # photo's left, bottom, width (fractions of the map; height follows the photo)
# ------------------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))  # files are found next to this script, from any working directory
UTM = 32715  # UTM zone 15S, metres

land = gpd.read_file(f"{HERE}/galapagos_land.geojson")
land["geometry"] = [Polygon(g.exterior) for g in land.geometry]  # drop lake holes


def coastal_route():
    """Española -> SE Santa Cruz -> along the east coast -> Daphne Major."""
    sc = land.to_crs(UTM).geometry[land.contains(Point(-90.36, -0.63))].iloc[0]
    ring = sc.buffer(COAST_GAP_KM * 1000).simplify(300).exterior
    ring = gpd.GeoSeries([ring], crs=UTM).to_crs(4326).iloc[0]
    xy = np.array(ring.coords)[:-1]
    if not ring.is_ccw:
        xy = xy[::-1]
    i = np.argmin(np.hypot(*(xy - ROUTE_ENTRY).T))
    j = np.argmin(np.hypot(*(xy - ROUTE_EXIT).T))
    coast = xy[i:j + 1] if i <= j else np.vstack([xy[i:], xy[:j + 1]])
    pts = np.vstack([ESPANOLA, coast, DAPHNE])
    for _ in range(3):                       # Chaikin corner smoothing
        q = 0.75 * pts[:-1] + 0.25 * pts[1:]
        r = 0.25 * pts[:-1] + 0.75 * pts[1:]
        pts = np.vstack([pts[0], np.column_stack([q, r]).reshape(-1, 2), pts[-1]])
    return pts, coast[0]


route, entry = coastal_route()
GEOD = Geod(ellps="WGS84")
ROUTE_KM = GEOD.line_length(route[:, 0], route[:, 1]) / 1000        # along the red dashed line
STRAIGHT_KM = GEOD.inv(*ESPANOLA, *DAPHNE)[2] / 1000              # Española - Daphne Major in a straight line


def draw_route(ax, lw, gap_deg=0.022):
    # stop the line just short of Daphne Major so the arrowhead does not cover it
    line = LineString(route)
    short = line.length - gap_deg
    body = np.array(LineString([line.interpolate(d).coords[0]
                                for d in np.linspace(0, short, 400)]).coords)
    ax.plot(body[:-5, 0], body[:-5, 1], color=RED, lw=lw, zorder=6,
            linestyle=(0, ROUTE_DASH), dash_capstyle="butt")
    ax.annotate("", xy=body[-1], xytext=body[-6], zorder=6,
                arrowprops=dict(arrowstyle="-|>,head_length=%g,head_width=%g" % ARROW_HEAD,
                                color=RED, lw=lw, shrinkA=0, shrinkB=0, mutation_scale=8 * lw / ROUTE_LW))


def draw_map(ax, scale=1.0):
    """The whole map on one axis. scale multiplies every text, line and marker size: 1 for a map FIG_WIDTH inches wide, larger for a
    larger map (the map's width in inches / FIG_WIDTH keeps the same look)."""
    k = scale
    halo = [pe.withStroke(linewidth=2 * k, foreground="white", alpha=0.55)]
    x0, x1, y0, y1 = EXTENT
    ax.set_facecolor(OCEAN)
    land.plot(ax=ax, facecolor=LAND, edgecolor=COAST, linewidth=0.5 * k, zorder=2)
    ax.set_aspect("equal"); ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.grid(False)

    # island labels: name, lon, lat, horizontal alignment, vertical alignment
    labels = [
        ("Isabela",       -91.10, -0.78, "center", "center"),
        ("Fernandina",    -91.52, -0.39, "center", "center"),
        ("Santiago",      -90.71, -0.26, "center", "center"),
        ("Santa Cruz",    -90.358, -0.627, "center", "center"),   # the middle of the island
        ("San Cristóbal", -89.430, -0.827, "center", "center"),   # the middle of the island
        ("Floreana",      -90.43, -1.29, "center", "center"),
        ("Española",      -89.74, -1.42, "right",  "top"),
        ("Santa Fe",      -89.99, -0.82, "left",   "center"),
        ("Baltra",        -90.21, -0.43, "left",   "center"),
        ("Pinzón",        -90.67, -0.67, "center", "top"),
        ("Rábida",        -90.71, -0.45, "center", "top"),
    ]
    for name, lon, lat, ha, va in labels:
        ax.text(lon, lat, name, ha=ha, va=va, fontsize=8 * k, color=TEXT,
                path_effects=halo, zorder=5)

    # Daphne Major: too small to see at this scale, so mark it
    ax.plot(*DAPHNE, "o", ms=6 * k, mfc="none", mec="black", mew=1.0 * k, zorder=7)
    ax.text(DAPHNE[0] - 0.05, DAPHNE[1] + 0.06, "Daphne Major", ha="left",
            va="bottom", fontsize=8.5 * k, fontweight="bold", color=TEXT,
            path_effects=halo, zorder=6)

    # dispersal route
    ax.plot(*ESPANOLA, "o", ms=3.5 * k, color=RED, zorder=6)
    draw_route(ax, ROUTE_LW * k)
    mid = (np.array(ESPANOLA) + entry) / 2
    ang = np.degrees(np.arctan2(entry[1] - ESPANOLA[1], ESPANOLA[0] - entry[0]))
    ax.text(mid[0] + 0.012, mid[1] + 0.012, DIST_LABEL or f"{ROUTE_KM:.0f} km", color=RED, fontsize=8 * k,
            rotation=-ang, rotation_mode="anchor", ha="center", va="bottom", zorder=6)

    if SHOW_TOWNS:
        towns = [("Puerto Ayora", -90.313, -0.743, "right", "top"),
                 ("Puerto Baquerizo Moreno", -89.610, -0.902, "right", "top"),
                 ("Puerto Villamil", -90.966, -0.956, "left", "top"),
                 ("Puerto Velasco Ibarra", -90.488, -1.275, "right", "bottom")]
        for name, lon, lat, ha, va in towns:
            ax.plot(lon, lat, "s", ms=2.5 * k, color=TEXT, zorder=6)
            ax.text(lon + (0.02 if ha == "left" else -0.02), lat, name, ha=ha, va=va,
                    fontsize=6.5 * k, style="italic", color=TEXT,
                    path_effects=halo, zorder=6)

    # scale bar (1 degree of longitude at the equator = 111.32 km)
    sb_km, sx, sy = 50, x0 + 0.08, y0 + 0.08
    ax.plot([sx, sx + sb_km / 111.32], [sy, sy], color="black", lw=1.5 * k,
            solid_capstyle="butt", zorder=6)
    ax.text(sx + sb_km / 111.32 / 2, sy + 0.025, f"{sb_km} km", ha="center",
            va="bottom", fontsize=7.5 * k, color="black", zorder=6)

    # north arrow
    ax.annotate("N", xy=(x0 + 0.12, y1 - 0.06), xytext=(x0 + 0.12, y1 - 0.28),
                ha="center", va="center", fontsize=8 * k, color="black",
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1 * k, mutation_scale=8 * k), zorder=6)

    ax.text(-91.45, -1.27, "Pacific Ocean", fontsize=10 * k, style="italic",
            color=SEA_TEXT, ha="center")
    ax.text(x1 - 0.09, y1 - 0.09, TITLE, fontsize=12 * k, style="italic", color=SEA_TEXT, ha="right", va="top")

    # axes
    ax.set_xticks([-91.5, -91, -90.5, -90, -89.5])
    ax.set_xticklabels(["91.5°W", "91°W", "90.5°W", "90°W", "89.5°W"])
    ax.set_yticks([-1, -0.5, 0])
    ax.set_yticklabels(["1°S", "0.5°S", "0°"])
    ax.tick_params(bottom=True, left=True, length=2.5 * k, width=0.5 * k, labelsize=7 * k, colors="#444444")
    for s in ax.spines.values():
        s.set_visible(True); s.set_linewidth(0.6 * k); s.set_color("#444444")

    # photo of the founder's species in the open ocean at the bottom right, labeled below
    photo = plt.imread(f"{HERE}/{BIRD_IMAGE}")
    bx, by, bw = BIRD_BOX
    bh = bw * photo.shape[0] / photo.shape[1] * (x1 - x0) / (y1 - y0)  # keep the photo's aspect ratio (the map's axes are x1-x0 by y1-y0)
    ax_bird = ax.inset_axes([bx, by, bw, bh])
    ax_bird.imshow(photo, interpolation="lanczos")
    ax_bird.set_xticks([]); ax_bird.set_yticks([]); ax_bird.grid(False)
    for s in ax_bird.spines.values():
        s.set_visible(True); s.set_linewidth(0.6 * k); s.set_color("#444444")
    for text, style, below in [("G. conirostris", "italic", 2.5), ("(founder)", "normal", 12.5)]:  # two lines, each this many points below the photo
        ax.annotate(text, (bx + bw / 2, by), xycoords="axes fraction", xytext=(0, -below * k), textcoords="offset points", style=style,
                    ha="center", va="top", fontsize=7.5 * k, color=TEXT, zorder=6)


if __name__ == "__main__":
    print(f"Straight-line distance: {STRAIGHT_KM:.1f} km")
    print(f"Length of drawn route:  {ROUTE_KM:.1f} km")
    set_theme()
    x0, x1, y0, y1 = EXTENT
    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_WIDTH * (y1 - y0) / (x1 - x0)))
    draw_map(ax)
    fig.savefig(f"{HERE}/galapagos_map.svg", dpi=600, bbox_inches="tight")
