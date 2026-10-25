#!/usr/bin/env python3
"""Figure 1 sandbox (a: MDFS schematic, b: all directed pairs). Writes only next to this script.

usage: python make_figure1.py [dir with {taxa,functions,combined}_fold*.tsv.gz] [outdir]
default data dir: ../figure1b_all_pairs/directed
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "figure1b_all_pairs", "directed")
OUT = sys.argv[2] if len(sys.argv) > 2 else HERE

INK, MUTED, GRID = "#1B2330", "#6A707A", "#E8EBEF"
AMBER = ("#FBE3A8", "#D99A1E", "#6B4500")
BLUE = ("#CFE0FA", "#4C7FD1", "#1E3F7A")
GREY = ("#EEF0F3", "#C5C9D0", "#6A707A")
GREEN = ("#DDF1E8", "#3E9B77", "#1F5E47")
NEVER = "#C9CDD3"
cmap = LinearSegmentedColormap.from_list("syn", ["#E1F2D8", "#86CFA8", "#2A9D8F", "#1B6A86"])
plt.rcParams.update({"font.family": ["Helvetica", "Arial", "DejaVu Sans"], "svg.fonttype": "none"})


def load(dt):
    parts = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(DATA, f"{dt}_fold*.tsv.gz")))]
    a = pd.concat(parts).groupby(["base", "contributing"], as_index=False).sum(numeric_only=True)
    a["x"] = a.base_1D / a.n
    a["y"] = a.total / a.n
    a["pct_syn"] = 100 * a.n_syn / a.n
    return a


# ---------- panel a (drawn in a 680 x 330 coordinate box, y down) ----------
def box(ax, x, y, w, h, colors, text=None, bold=False, r=8, fs=9.5, lw=0.9, txt_dy=0):
    fill, edge, ink = colors
    ax.add_patch(FancyBboxPatch((x + r / 2, y + r / 2), w - r, h - r, boxstyle=f"round,pad={r / 2},rounding_size={r}",
                                fc=fill, ec=edge, lw=lw))
    if text:
        ax.text(x + w / 2, y + h / 2 + txt_dy, text, ha="center", va="center", fontsize=fs,
                fontweight="bold" if bold else "normal", color=ink)


def arrow(ax, p, q, strong=True):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9, shrinkA=0, shrinkB=0,
                                 color=INK if strong else "#B4BAC4", lw=1.6 if strong else 0.9,
                                 ls="-" if strong else (0, (3, 3))))


def outcome(ax, x, y, w, h):
    box(ax, x, y, w, h, ("#FFFFFF", INK, INK), r=8, lw=1.2)
    ax.text(x + w / 2, y + h / 2 - 12, "Outcome", ha="center", va="center", fontsize=9.5, fontweight="bold", color=INK)
    ax.text(x + w / 2, y + h / 2 + 5, "disease", ha="center", va="center", fontsize=8.5, color=MUTED)
    ax.text(x + w / 2, y + h / 2 + 17, "label", ha="center", va="center", fontsize=8.5, color=MUTED)


def panel_a(ax):
    ax.set_xlim(14, 684); ax.set_ylim(322, 10); ax.axis("off")
    for x0, w, title, sub in ((30, 300, "1D mode", "each feature tested alone"),
                              (350, 320, "2D mode", "feature pairs tested jointly")):
        ax.add_patch(FancyBboxPatch((x0 + 6, 18), w - 12, 294, boxstyle="round,pad=6,rounding_size=10",
                                    fc="#F7F8FA", ec="#E1E5EA", lw=0.8))
        ax.text(x0 + 16, 38, title, fontsize=11.5, fontweight="bold", color=INK, va="center")
        ax.text(x0 + 16, 56, sub, fontsize=8.5, color=MUTED, va="center")
    feats = [("Taxon 1", AMBER, True), ("Taxon 2", GREY, False), ("Taxon 3", GREY, False),
             ("Function 2", GREY, False), ("Function 1", BLUE, True)]
    for i, (name, col, sig) in enumerate(feats):
        y = 72 + 44 * i
        box(ax, 46, y, 124, 34, col, name, bold=sig)
        ty = (152, 162, 175, 188, 198)[i]
        arrow(ax, (172, y + 17), (238, ty), strong=sig)
    outcome(ax, 240, 140, 76, 70)
    for gy, pair in ((72, (("Taxon 1", AMBER, True), ("Taxon 2", GREY, False))),
                     (196, (("Taxon 3", GREY, False), ("Function 2", GREY, False)))):
        box(ax, 366, gy, 146, 108, GREEN, r=12)
        for j, (name, col, sig) in enumerate(pair):
            box(ax, 378, gy + 10 + 38 * j, 122, 30, col, name, bold=sig)
        ax.text(439, gy + 95, "synergistic pair", ha="center", va="center", fontsize=8.5, fontweight="bold",
                color=GREEN[2])
    arrow(ax, (514, 126), (574, 170))
    arrow(ax, (514, 250), (574, 204))
    outcome(ax, 576, 152, 82, 70)


# ---------- panel b ----------
def panel_b(axes, cax):
    data = {dt: load(dt) for dt in ("taxa", "functions", "combined")}
    lim = max(max(d.x.max(), d.y.max()) for d in data.values()) * 1.03
    for ax, (dt, d) in zip(axes, data.items()):
        never = d[d.n_syn == 0]
        syn = d[d.n_syn > 0].sort_values("pct_syn")      # most consistent drawn on top
        ax.scatter(never.x, never.y, s=3, c=NEVER, alpha=0.45, lw=0, rasterized=True)
        sc = ax.scatter(syn.x, syn.y, s=2.5, c=syn.pct_syn, cmap=cmap, vmin=0, vmax=100, alpha=0.5, lw=0,
                        rasterized=True)
        ax.plot([0, lim], [0, lim], ls="--", lw=0.8, c=INK)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
        ax.set_title(f"{dt.capitalize()}", loc="left", fontsize=10.5, fontweight="bold", color=INK, pad=14)
        ax.text(0, 1.02, f"{len(d):,} pairs  |  {len(syn):,} synergistic in ≥1 run", transform=ax.transAxes,
                fontsize=7.5, color=MUTED)
        ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color("#9AA1AB")
        ax.tick_params(colors=MUTED, labelsize=8, length=3)
        ax.set_xlabel("Base feature 1D IG", fontsize=9, color=MUTED)
    axes[0].set_ylabel("Base 1D IG + 2D IG added by contributing feature", fontsize=9, color=MUTED)
    cb = plt.colorbar(sc, cax=cax)
    cb.set_label("% of fold-runs in which pair is synergistic", fontsize=8.5, color=MUTED)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8, colors=MUTED, length=2)
    axes[2].plot([], [], ls="--", lw=0.8, c=INK, label="y = x (no 2D gain)")
    axes[2].scatter([], [], s=12, c=NEVER, label="never synergistic")
    axes[2].legend(loc="lower right", fontsize=7.5, frameon=False, labelcolor=MUTED)


fig = plt.figure(figsize=(10.5, 7.0), facecolor="white")
gs = fig.add_gridspec(2, 3, height_ratios=[0.8, 1.3], hspace=0.02, wspace=0.14,
                      left=0.065, right=0.925, top=0.97, bottom=0.075)
axa = fig.add_subplot(gs[0, :])
axa.set_aspect("auto")
panel_a(axa)
axb = [fig.add_subplot(gs[1, i]) for i in range(3)]
cax = fig.add_axes([0.945, 0.20, 0.012, 0.28])
panel_b(axb, cax)
fig.text(0.005, 0.975, "a", fontsize=16, fontweight="bold", color=INK, va="top")
fig.text(0.005, axb[0].get_position().y1 + 0.05, "b", fontsize=16, fontweight="bold", color=INK, va="top")
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"figure1_combined.{ext}"), dpi=200)
print("saved")
