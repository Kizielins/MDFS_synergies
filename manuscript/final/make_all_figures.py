#!/usr/bin/env python3
"""Main-text Figures 1-4 of the MDFS synergy manuscript, one script, one shared style.

Built from the existing plotting code:
  Figure 1  figure1_sandbox.py / R1_synergistic_pairs/figure1_combined/make_figure1.py   (a: schematic, b: all directed pairs)
  Figure 2  make_figure2.py / R2_classification_performance/effect_sizes.ipynb           (hub pairs: selection frequency + effect sizes)
  Figure 3  R3_pair_transformations/transformation_preference.ipynb                       (LR vs GM, 3 x 3 pairs)
  Figure 4  R4_broader_application/analyse_synergy_gains.ipynb                            (non-CRC cohorts, 4 panels)

Shared style (identical for all figures)
  * one sans-serif family (Arial > Helvetica > Liberation Sans > DejaVu Sans), also used for the italic/bold
    parts of species names, so every glyph in every figure comes from the same font
  * panel labels: lower-case, bold, 8 pt, no brackets (a, b, c ...), always placed at the top-left corner of the
    panel they belong to
  * text sizes (pt, at the final print size): panel label 8, axis title / panel title 7, tick labels 6,
    annotations / legends 6 (nothing below 5.5 pt)
  * figure width 180 mm (two columns). Heights are fixed per figure, nothing is cropped on saving
    (no bbox_inches="tight"), so the files have exactly the size they will be printed at
  * no title / "Figure N" label inside the image (legends go into the manuscript text)
  * output: vector PDF (fonts embedded as TrueType, dense scatter layers rasterised at --dpi) and PNG at --dpi
    (default 600 dpi, i.e. 4252 px across); --formats pdf,png,tiff selects the formats

Usage (all paths are examples; every figure is optional, use --only to run a subset)
  python make_all_figures.py --out high_res \
      --fig1-data  <dir with {taxa,functions,combined}_fold*.tsv.gz> \
      --fig2-data  <repo>/R2_classification_performance \
      --prepared-data <repo>/M1_data_preparation/prepared_data  --fig3-candidates <repo>/R3_pair_transformations/figure4_candidates.tsv \
      --fig4-data  <repo>/R4_broader_application/synergy_by_cohort.tsv
  python make_all_figures.py --only 2,4 --fig2-data ... --fig4-data ...
Figure 3 needs the prepared cohort matrices once; the per-cohort effect sizes it computes are cached in
--fig3-cache (default: <out>/figure3_effects.csv) so later runs do not need them.
"""
import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch  # noqa: E402

# ----------------------------------------------------------------------------------------------- shared style
FONT_STACK = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
INK, MUTED, GRID, AXIS = "#1B2330", "#6A707A", "#E8EBEF", "#9AA1AB"
C_SIG, C_NS = "#1B6A86", "#C5C9D0"

W_IN = 180 / 25.4            # figure width: 180 mm
FS_PANEL, FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 8, 7, 7, 6, 6

plt.rcParams.update({
    "font.family": FONT_STACK, "font.size": FS_TICK,
    "mathtext.fontset": "custom", "mathtext.rm": FONT_STACK[0], "mathtext.it": FONT_STACK[0] + ":italic",
    "mathtext.bf": FONT_STACK[0] + ":bold", "mathtext.default": "regular", "mathtext.cal": FONT_STACK[0], "mathtext.sf": FONT_STACK[0], "mathtext.tt": FONT_STACK[0],
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.unicode_minus": True,
})
# make the mathtext font fall back to the first installed family of the stack
from matplotlib import font_manager as _fm  # noqa: E402
_installed = {f.name for f in _fm.fontManager.ttflist}
_fam = next((f for f in FONT_STACK if f in _installed), "DejaVu Sans")
plt.rcParams.update({"font.family": [_fam], "mathtext.rm": _fam, "mathtext.it": _fam + ":italic",
                     "mathtext.bf": _fam + ":bold", "mathtext.cal": _fam, "mathtext.sf": _fam, "mathtext.tt": _fam})


def style_axes(ax, grid="both"):
    """Axis look shared by every plotted panel."""
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelsize=FS_TICK)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.5)
    ax.set_axisbelow(True)


def axis_label(ax, x=None, y=None):
    if x is not None:
        ax.set_xlabel(x, fontsize=FS_LABEL, color=MUTED)
    if y is not None:
        ax.set_ylabel(y, fontsize=FS_LABEL, color=MUTED)


def panel_label(fig, letter, x, y):
    """Lower-case bold panel letter at figure coordinates (x, y) = top-left of the letter."""
    fig.text(x, y, letter, fontsize=FS_PANEL, fontweight="bold", color=INK, ha="left", va="top")


def label_at_axes(fig, ax, letter, dx_pt=0.0, dy_pt=2.0):
    """Panel letter just above the top-left corner of the tight bounding box (ticks, axis title, panel title
    included) of `ax`, so that it can never collide with tick labels or axis titles."""
    fig.canvas.draw()
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bb.x0 + dx_pt / (fig.get_figwidth() * 72), bb.y1 + dy_pt / (fig.get_figheight() * 72), letter,
             fontsize=FS_PANEL, fontweight="bold", color=INK, ha="left", va="bottom")


def labels_row(fig, items, dy_pt=3.0):
    """Panel letters for several panels that share a row: all letters sit at the same height (just above the highest
    panel extent of the row) and at the left edge of their own panel extent."""
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    boxes = [(ax.get_tightbbox(fig.canvas.get_renderer()).transformed(inv), l) for ax, l in items]
    ytop = max(b.y1 for b, _ in boxes) + dy_pt / (fig.get_figheight() * 72)
    for b, l in boxes:
        fig.text(b.x0, ytop, l, fontsize=FS_PANEL, fontweight="bold", color=INK, ha="left", va="bottom")


def species(name):
    """'Parvimonas_micra' -> mathtext italic 'P. micra'; other feature names (pathways) are returned unchanged."""
    s = str(name)
    if ":" in s or "-PWY" in s or s.upper() == s:
        return s.split(":", 1)[0].replace("_", " ")
    parts = s.replace("_", " ").split(" ")
    if len(parts) >= 2:
        rest = " ".join(parts[1:]).replace(" ", "\\ ")
        return "$\\mathit{%s.\\ %s}$" % (parts[0][0], rest)
    return s


def save(fig, out, name, formats, dpi):
    os.makedirs(out, exist_ok=True)
    for ext in formats:
        path = os.path.join(out, f"{name}.{ext}")
        kw = {"dpi": dpi, "facecolor": "white"}
        if ext == "tiff":
            kw["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(path, **kw)
        print("  wrote", path)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------- Figure 1
AMBER = ("#FBE3A8", "#D99A1E", "#6B4500")
BLUE = ("#CFE0FA", "#4C7FD1", "#1E3F7A")
GREY = ("#EEF0F3", "#C5C9D0", "#6A707A")
GREEN = ("#DDF1E8", "#3E9B77", "#1F5E47")
NEVER = "#C9CDD3"
CMAP1 = LinearSegmentedColormap.from_list("syn", ["#E1F2D8", "#86CFA8", "#2A9D8F", "#1B6A86"])


def _load_pairs(data, dt):
    parts = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(data, f"{dt}_fold*.tsv.gz")))]
    if not parts:
        sys.exit(f"Figure 1: no {dt}_fold*.tsv.gz files in {data}")
    a = pd.concat(parts).groupby(["base", "contributing"], as_index=False).sum(numeric_only=True)
    a["x"] = a.base_1D / a.n
    a["y"] = a.total / a.n
    a["pct_syn"] = 100 * a.n_syn / a.n
    return a


def figure1(data, out, formats, dpi):
    # ---- panel a: schematic, drawn in a 680 x 330 box (y down); text sizes 6 pt (box labels) / 5.5 pt (sub text)
    def box(ax, x, y, w, h, colors, text=None, bold=False, r=8, lw=0.6):
        fill, edge, ink = colors
        ax.add_patch(FancyBboxPatch((x + r / 2, y + r / 2), w - r, h - r, boxstyle=f"round,pad={r / 2},rounding_size={r}",
                                    fc=fill, ec=edge, lw=lw))
        if text:
            ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=FS_TICK + 0.5,
                    fontweight="bold" if bold else "normal", color=ink)

    def arrow(ax, p, q, strong=True):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=6, shrinkA=0, shrinkB=0,
                                     color=INK if strong else "#B4BAC4", lw=1.0 if strong else 0.6,
                                     ls="-" if strong else (0, (3, 3))))

    def outcome(ax, x, y, w, h):
        box(ax, x, y, w, h, ("#FFFFFF", INK, INK), r=8, lw=0.8)
        ax.text(x + w / 2, y + h / 2 - 14, "Outcome", ha="center", va="center", fontsize=FS_TICK + 0.5,
                fontweight="bold", color=INK)
        ax.text(x + w / 2, y + h / 2 + 6, "disease", ha="center", va="center", fontsize=FS_ANNOT, color=MUTED)
        ax.text(x + w / 2, y + h / 2 + 20, "label", ha="center", va="center", fontsize=FS_ANNOT, color=MUTED)

    def panel_a(ax):
        ax.set_xlim(14, 684); ax.set_ylim(322, 10); ax.axis("off")
        for x0, w, title, sub in ((30, 300, "1D mode", "each feature tested alone"),
                                  (350, 320, "2D mode", "feature pairs tested jointly")):
            ax.add_patch(FancyBboxPatch((x0 + 6, 18), w - 12, 294, boxstyle="round,pad=6,rounding_size=10",
                                        fc="#F7F8FA", ec="#E1E5EA", lw=0.6))
            ax.text(x0 + 16, 38, title, fontsize=FS_TITLE, fontweight="bold", color=INK, va="center")
            ax.text(x0 + 16, 56, sub, fontsize=FS_ANNOT, color=MUTED, va="center")
        feats = [("Taxon 1", AMBER, True), ("Taxon 2", GREY, False), ("Taxon 3", GREY, False),
                 ("Function 2", GREY, False), ("Function 1", BLUE, True)]
        for i, (name, col, sig) in enumerate(feats):
            y = 72 + 44 * i
            box(ax, 46, y, 124, 34, col, name, bold=sig)
            arrow(ax, (172, y + 17), (238, (152, 162, 175, 188, 198)[i]), strong=sig)
        outcome(ax, 240, 140, 76, 70)
        for gy, pair in ((72, (("Taxon 1", AMBER, True), ("Taxon 2", GREY, False))),
                         (196, (("Taxon 3", GREY, False), ("Function 2", GREY, False)))):
            box(ax, 366, gy, 146, 108, GREEN, r=12)
            for j, (name, col, sig) in enumerate(pair):
                box(ax, 378, gy + 10 + 38 * j, 122, 30, col, name, bold=sig)
            ax.text(439, gy + 95, "synergistic pair", ha="center", va="center", fontsize=FS_ANNOT, fontweight="bold",
                    color=GREEN[2])
        arrow(ax, (514, 126), (574, 170))
        arrow(ax, (514, 250), (574, 204))
        outcome(ax, 576, 152, 82, 70)

    fig = plt.figure(figsize=(W_IN, 4.6), facecolor="white")
    gs = fig.add_gridspec(2, 3, height_ratios=[0.8, 1.3], hspace=0.10, wspace=0.28,
                          left=0.075, right=0.885, top=0.955, bottom=0.085)
    axa = fig.add_subplot(gs[0, :])
    panel_a(axa)
    axb = [fig.add_subplot(gs[1, i]) for i in range(3)]
    cax = fig.add_axes([0.912, 0.20, 0.012, 0.28])

    d = {dt: _load_pairs(data, dt) for dt in ("taxa", "functions", "combined")}
    lim = max(max(x.x.max(), x.y.max()) for x in d.values()) * 1.03
    sc = None
    for ax, (dt, df) in zip(axb, d.items()):
        never = df[df.n_syn == 0]
        syn = df[df.n_syn > 0].sort_values("pct_syn")
        ax.scatter(never.x, never.y, s=1.2, c=NEVER, alpha=0.45, lw=0, rasterized=True)
        sc = ax.scatter(syn.x, syn.y, s=1.0, c=syn.pct_syn, cmap=CMAP1, vmin=0, vmax=100, alpha=0.5, lw=0, rasterized=True)
        ax.plot([0, lim], [0, lim], ls="--", lw=0.6, c=INK)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
        ax.set_title(dt.capitalize(), loc="left", fontsize=FS_TITLE, fontweight="bold", color=INK, pad=11)
        ax.text(0, 1.02, f"{len(df):,} pairs  |  {len(syn):,} synergistic in ≥1 run", transform=ax.transAxes,
                fontsize=FS_ANNOT - 0.5, color=MUTED)
        style_axes(ax)
        axis_label(ax, x="Base feature 1D IG")
    axis_label(axb[0], y="Base 1D IG + 2D IG added by contributing feature")
    cb = plt.colorbar(sc, cax=cax)
    cb.set_label("% of fold-runs in which pair is synergistic", fontsize=FS_LABEL, color=MUTED)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=FS_TICK, colors=MUTED, length=2)
    axb[2].plot([], [], ls="--", lw=0.6, c=INK, label="y = x (no 2D gain)")
    axb[2].scatter([], [], s=8, c=NEVER, label="never synergistic")
    axb[2].legend(loc="lower right", fontsize=FS_ANNOT, frameon=False, labelcolor=MUTED)
    panel_label(fig, "a", 0.005, 0.995)
    label_at_axes(fig, axb[0], "b", dy_pt=6)
    save(fig, out, "Figure_1", formats, dpi)


# ---------------------------------------------------------------------------------------------------- Figure 2
def figure2(data, out, formats, dpi):
    C_BASE, C_CONTRIB, C_PAIR = "#B9C7DD", "#F2C777", "#3E9B77"
    tests = pd.read_csv(os.path.join(data, "cohort_level_tests.csv"))
    sel = pd.read_csv(os.path.join(data, "final_synthetic_features_table.csv"))
    meta = pd.read_csv(os.path.join(data, "meta_results_all_pairs.csv"))
    ranked = meta.merge(sel[["Data Type", "Feature", "Selection % Gain"]].rename(
        columns={"Data Type": "data_type", "Feature": "feature_pair"}), on=["data_type", "feature_pair"], how="left")
    pm = ranked[ranked.feature_pair.str.contains("Parvimonas_micra")]
    top = pm[pm["Selection % Gain"] >= 99.0]
    pool = top if len(top) else pm
    r1 = pool.loc[pool.z_meta.abs().idxmax()]
    fn = ranked[ranked.feature_pair.str.contains("Fusobacterium_nucleatum") & (ranked.feature_pair != r1.feature_pair)]
    r2 = fn.loc[fn.z_meta.abs().idxmax()]
    pairs = [(r1.data_type, r1.feature_pair), (r2.data_type, r2.feature_pair)]

    def cohort_label(c):
        c = c.replace("Public_study__", "").replace("This_study_", "")
        if c.startswith("cohort"):
            c = c.split("__", 1)[1]
        return c.replace("_", " ")

    def sp_full(n):
        g, s = n.replace("_", " ").split(" ", 1)
        return f"$\\mathit{{{g}}}$\n$\\mathit{{{s}}}$"

    fig, axes = plt.subplots(2, 2, figsize=(W_IN, 5.6), gridspec_kw={
        "width_ratios": [1.2, 1.1], "hspace": 0.85, "wspace": 0.62, "left": 0.115, "right": 0.985, "top": 0.885, "bottom": 0.075})
    for row, (dt, fp) in enumerate(pairs):
        pd_ = tests[(tests.data_type == dt) & (tests.feature_pair == fp)].copy()
        s = sel[(sel["Data Type"] == dt) & (sel.Feature == fp)].iloc[0]
        f = sorted([(s.Feature1, s["Feature1_Selection%"]), (s.Feature2, s["Feature2_Selection%"])], key=lambda t: t[1])
        (bn, bv), (cn, cv) = f
        pv = s["Selection %"]
        avg = np.mean([bv, cv])

        ax = axes[row, 0]
        ax.bar([0, 1, 2], [bv, cv, pv], width=0.62, color=[C_BASE, C_CONTRIB, C_PAIR], zorder=3)
        ax.axhline(avg, color=INK, ls="--", lw=0.7, zorder=4, label=f"mean of single features ({avg:.0f}%)")
        ax.legend(loc="upper left", frameon=False, fontsize=FS_ANNOT, labelcolor=MUTED, handlelength=1.6)
        for x, v in enumerate([bv, cv, pv]):
            ax.text(x, v + 1.5, f"{v:.0f}%", ha="center", va="bottom", fontsize=FS_ANNOT, color=INK, fontweight="bold")
        ax.annotate("", xy=(2.43, pv), xytext=(2.43, avg), arrowprops=dict(arrowstyle="<->", lw=0.7, color=INK))
        ax.text(2.5, (pv + avg) / 2, f"+{pv - avg:.0f} pp", ha="left", va="center", fontsize=FS_ANNOT, fontweight="bold",
                color=C_PAIR)
        ax.set_xlim(-0.5, 3.1); ax.set_ylim(0, 125); ax.set_yticks(range(0, 101, 20)); ax.set_xticks([0, 1, 2])
        ax.set_xticklabels([sp_full(bn) + "\n(base)", sp_full(cn) + "\n(contributor)", "Synthetic\npair"], fontsize=FS_TICK - 0.5)
        axis_label(ax, y="Selection frequency (% of LOCO folds)")
        ax.set_title("Selection frequency", loc="left", fontsize=FS_TITLE, fontweight="bold", color=INK, pad=6)
        style_axes(ax, grid="y")

        ax2 = axes[row, 1]
        pd_ = pd_.assign(name=pd_.cohort.map(cohort_label)).sort_values("rbc", ascending=True)
        y = np.arange(len(pd_))
        sig = pd_.p_value_fdr.fillna(1) < 0.05
        ax2.barh(y, pd_.rbc, height=0.62, color=np.where(sig, C_SIG, C_NS), zorder=3)
        ax2.axvline(0, color=INK, lw=0.7, zorder=4)
        ax2.set_yticks(y); ax2.set_yticklabels(pd_.name, fontsize=FS_TICK)
        axis_label(ax2, x="Rank-biserial r (CRC vs control)")
        n_sig = int(sig.sum()); n_cons = int((pd_.rbc > 0).sum() if pd_.rbc.median() > 0 else (pd_.rbc < 0).sum())
        ax2.set_title("Effect size per cohort", loc="left", fontsize=FS_TITLE, fontweight="bold", color=INK, pad=14)
        ax2.text(0, 1.03, f"{n_sig}/{len(pd_)} cohorts significant (FDR q < 0.05), {n_cons}/{len(pd_)} same direction",
                 transform=ax2.transAxes, fontsize=FS_ANNOT, color=MUTED)
        style_axes(ax2, grid="x")
        ax2.legend(handles=[Patch(fc=C_SIG, label="FDR q < 0.05"), Patch(fc=C_NS, label="FDR q ≥ 0.05")],
                   loc="lower right", frameon=False, fontsize=FS_ANNOT, labelcolor=MUTED, handlelength=1.2)
        ax.text(0.0, 1.30, f"Pair {row + 1}:  {species(bn)} (base) + {species(cn)} (contributor)", transform=ax.transAxes,
                fontsize=FS_TITLE, fontweight="bold", color=INK, va="top")
        print("  pair", row + 1, dt, fp, f"base={bv:.0f} contrib={cv:.0f} pair={pv:.0f} sig={n_sig}/{len(pd_)}")
    labels_row(fig, [(axes[0, 0], "a"), (axes[0, 1], "b")], dy_pt=4)
    save(fig, out, "Figure_2", formats, dpi)


# ---------------------------------------------------------------------------------------------------- Figure 3
def _bh(p):
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    po = p[ok]; o = np.argsort(po); m = po.size
    qr = np.minimum.accumulate((po[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    qo = np.empty_like(po); qo[o] = np.clip(qr, 0, 1); q[ok] = qo
    return q


def _find_col(name, cols):
    if name in cols:
        return name
    if isinstance(name, str) and ":" in name:
        pid = name.split(":")[0].strip()
        if pid in cols:
            return pid
        for c in cols:
            if str(c).startswith(pid):
                return c
    low = {str(c).lower(): c for c in cols}
    return low.get(str(name).lower())


def _select_fig3_pairs(cand_path):
    """Same automatic choice of the nine pairs as in transformation_preference.ipynb (read from
    figure4_candidates.tsv): 3 LR-exclusive, 3 GM-exclusive (taxa/combined), 3 functional pairs."""
    c = pd.read_csv(cand_path, sep="\t")
    c["pair_key"] = [tuple(sorted((a, b))) for a, b in zip(c.f1, c.f2)]
    c = c.sort_values("selection_pct", ascending=False, kind="stable").drop_duplicates("pair_key").sort_index()
    tx = c[c.data_type.isin(["taxa", "combined"])]
    r1 = tx[tx.types_ge_50pct == "LR-only"].sort_values("contrast_LR_minus_GM", ascending=False).head(3)
    r2 = tx[tx.types_ge_50pct == "GM-only"].sort_values("contrast_LR_minus_GM", ascending=True).head(3)
    fn = c[c.data_type == "functions"].copy()
    fn["pref_gm"] = (fn.types_ge_50pct == "GM-only").astype(int)
    r3 = fn.sort_values(["pref_gm", "median_abs_r_GM"], ascending=[False, False]).head(3)
    return pd.concat([r1, r2, r3])


def _compute_fig3_effects(sel, prepared):
    from scipy.stats import mannwhitneyu
    eps = 1e-9
    rows = []
    for r in sel.itertuples():
        stems = sorted(os.path.basename(p)[2:-len(f"_{r.data_type}.tsv")]
                       for p in glob.glob(os.path.join(prepared, f"X_*_{r.data_type}.tsv")))
        for ft in ("LR", "GM"):
            tmp = []
            for st in stems:
                X = pd.read_csv(os.path.join(prepared, f"X_{st}_{r.data_type}.tsv"), sep="\t", index_col=0)
                yp = os.path.join(prepared, f"y_{st}_{r.data_type}.tsv")
                if not os.path.exists(yp):
                    cands = sorted(glob.glob(os.path.join(prepared, f"y_{st}_*.tsv")))
                    if not cands:
                        continue
                    yp = cands[0]
                yd = pd.read_csv(yp, sep="\t", index_col=0)
                if "Study condition" not in yd.columns:
                    continue
                X = X.div(X.sum(axis=1), axis=0)
                y = yd["Study condition"].astype(int)
                idx = X.index.intersection(y.index); X, y = X.loc[idx], y.loc[idx]
                c1, c2 = _find_col(r.f1, X.columns), _find_col(r.f2, X.columns)
                if c1 is None or c2 is None:
                    continue
                a, b = X[c1].to_numpy(float), X[c2].to_numpy(float)
                v = np.log((a + eps) / (b + eps)) if ft == "LR" else np.sqrt((a + eps) * (b + eps))
                d, h = v[(y == 0).to_numpy()], v[(y == 1).to_numpy()]
                if len(d) == 0 or len(h) == 0:
                    continue
                U, p = mannwhitneyu(d, h, alternative="two-sided")
                tmp.append({"cohort": st.replace("Public_study__", "").replace("This_study_", ""),
                            "rbc": 2 * U / (len(d) * len(h)) - 1, "p": p})
            df = pd.DataFrame(tmp)
            if len(df):
                df["q"] = _bh(df.p.values)
                df.insert(0, "encoding", ft); df.insert(0, "feature_pair", r.feature_pair); df.insert(0, "data_type", r.data_type)
                rows.append(df)
    return pd.concat(rows)


def figure3(cand_path, prepared, cache, out, formats, dpi):
    sel = _select_fig3_pairs(cand_path)
    if os.path.exists(cache):
        eff = pd.read_csv(cache)
        print("  per-cohort effects read from", cache)
    else:
        if not prepared or not os.path.isdir(prepared):
            sys.exit("Figure 3: need --prepared-data (M1_data_preparation/prepared_data) or an existing --fig3-cache")
        eff = _compute_fig3_effects(sel, prepared)
        os.makedirs(os.path.dirname(os.path.abspath(cache)), exist_ok=True)
        eff.to_csv(cache, index=False)
        print("  per-cohort effects cached in", cache)

    LR_SIG, LR_NS, GM_SIG, GM_NS = "#2166ac", "#92c5de", "#b2182b", "#f4a582"
    fig = plt.figure(figsize=(W_IN, 6.9), facecolor="white")
    gs = fig.add_gridspec(3, 3, left=0.085, right=0.985, top=0.935, bottom=0.115, hspace=0.78, wspace=0.12)
    allr = eff.rbc.abs().max()
    xlim = min(1.0, max(0.25, float(allr) + 0.06))
    row_names = ["LR-exclusive", "GM-exclusive", "Function / strongest GM"]
    letters = "abcdefghi"
    first_axes = {}
    row_axes = {}
    for idx, r in enumerate(sel.itertuples()):
        rr, cc = divmod(idx, 3)
        ax = fig.add_subplot(gs[rr, cc])
        e = eff[(eff.data_type == r.data_type) & (eff.feature_pair == r.feature_pair)]
        lr = e[e.encoding == "LR"].set_index("cohort"); gm = e[e.encoding == "GM"].set_index("cohort")
        common = sorted(set(lr.index) & set(gm.index))
        lr, gm = lr.loc[common], gm.loc[common]
        dom = "LR" if np.median(np.abs(lr.rbc)) > np.median(np.abs(gm.rbc)) else "GM"
        order = np.argsort(-np.abs((lr if dom == "LR" else gm).rbc.values))
        lr, gm = lr.iloc[order], gm.iloc[order]
        y = np.arange(len(lr))
        for j in y:
            ax.plot([lr.rbc.iloc[j], gm.rbc.iloc[j]], [j, j], color="gray", lw=0.4, alpha=0.5, zorder=2)
        ax.scatter(lr.rbc, y, c=np.where(lr.q < 0.05, LR_SIG, LR_NS), s=14, marker="s", ec="black", lw=0.3, zorder=3)
        ax.scatter(gm.rbc, y, c=np.where(gm.q < 0.05, GM_SIG, GM_NS), s=14, marker="o", ec="black", lw=0.3, zorder=3)
        ax.axvline(0, color="black", ls="--", lw=0.7, alpha=0.75)
        ax.set_yticks([]); ax.invert_yaxis(); ax.set_xlim(-xlim, xlim)
        style_axes(ax, grid="x")
        ax.set_title(f"{r.data_type} | best: {dom} | n = {len(lr)}", fontsize=FS_ANNOT, fontweight="bold", color=INK,
                     pad=13, loc="center")
        ax.text(0.5, 1.02, f"{species(r.f1)} × {species(r.f2)}", transform=ax.transAxes, ha="center", va="bottom",
                fontsize=FS_ANNOT, color=INK)
        if cc == 0:
            ax.set_ylabel(row_names[rr], fontsize=FS_LABEL, color=MUTED, fontweight="bold", labelpad=4)
            first_axes[rr] = ax
        row_axes.setdefault(rr, []).append((ax, letters[idx]))
    for rr in sorted(row_axes):
        labels_row(fig, row_axes[rr], dy_pt=2)
    fig.supxlabel("Effect size r (CRC vs control)", fontsize=FS_LABEL, color=MUTED, y=0.065)
    handles = [Line2D([0], [0], marker="s", color="w", markerfacecolor=LR_SIG, markeredgecolor="black", markersize=5, label="LR (FDR q < 0.05)"),
               Line2D([0], [0], marker="s", color="w", markerfacecolor=LR_NS, markeredgecolor="black", markersize=5, label="LR (n.s.)"),
               Line2D([0], [0], marker="o", color="w", markerfacecolor=GM_SIG, markeredgecolor="black", markersize=5, label="GM (FDR q < 0.05)"),
               Line2D([0], [0], marker="o", color="w", markerfacecolor=GM_NS, markeredgecolor="black", markersize=5, label="GM (n.s.)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=FS_ANNOT, labelcolor=MUTED,
               bbox_to_anchor=(0.5, 0.0), columnspacing=2.0, handletextpad=0.3)
    save(fig, out, "Figure_3", formats, dpi)


# ---------------------------------------------------------------------------------------------------- Figure 4
# Two blue arrows that were drawn by hand on panel a of the submitted figure (start -> head, axes fractions).
FIG4_ARROWS = [((0.436, 0.12), (0.350, 0.47)), ((0.634, 0.29), (0.546, 0.64))]


def figure4(path, out, formats, dpi, arrows=True):
    import seaborn as sns
    from scipy.stats import gaussian_kde
    CRC = {"VogtmannE_2016", "WirbelJ_2018", "YachidaS_2019", "YuJ_2015", "ZellerG_2014"}
    s = pd.read_csv(path, sep="\t")
    s = s[~s.cohort.isin(CRC)].copy()
    s["pair_key"] = [" + ".join(sorted((str(a), str(b)))) for a, b in zip(s.base, s.contributing)]
    s = s.sort_values("total_IG_mean", ascending=False).drop_duplicates(["cohort", "pair_key"]).reset_index(drop=True)
    s["total_IG_mean"] = s.total_IG_mean.astype(float); s["base_IG_1D_mean"] = s.base_IG_1D_mean.astype(float)
    s = s.dropna(subset=["ig_gain_pct", "cohort", "total_IG_mean"])
    s["pct_med"] = s.total_IG_mean / s.groupby("cohort").total_IG_mean.transform("median") * 100
    s["log_gain"] = np.log10(s.ig_gain_pct.clip(lower=1))
    b_df, c_df = s[s.ig_gain_pct >= 50], s[s.total_IG_mean >= 5]
    cohorts = sorted(s.cohort.unique())
    pal = sns.color_palette("tab20", n_colors=len(cohorts))
    print("  non-CRC cohorts:", len(cohorts), "| unordered pairs:", len(s))

    fig, axes = plt.subplots(2, 2, figsize=(W_IN, 6.2), gridspec_kw={
        "hspace": 0.42, "wspace": 0.62, "left": 0.155, "right": 0.985, "top": 0.945, "bottom": 0.115})
    ax = axes[0, 0]
    x, y = s.base_IG_1D_mean, s.total_IG_mean
    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    sc = ax.scatter(x, y, c=s.log_gain, cmap="viridis", alpha=0.6, s=5, lw=0, rasterized=True)
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, label="x = y")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.legend(loc="upper left", fontsize=FS_ANNOT, frameon=True, framealpha=0.85, edgecolor="none", labelcolor=MUTED)
    axis_label(ax, x="Base feature 1D IG", y="Total IG (1D + 2D added)")
    style_axes(ax)
    cb = plt.colorbar(sc, ax=ax, pad=0.02, fraction=0.05)
    cb.set_label("log$_{10}$ IG gain %", fontsize=FS_LABEL, color=MUTED)
    cb.ax.tick_params(labelsize=FS_TICK, colors=MUTED, length=2); cb.outline.set_visible(False)
    if arrows:
        for (x0, y0), (x1, y1) in FIG4_ARROWS:
            ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=7,
                                         lw=0.8, color=C_SIG, shrinkA=0, shrinkB=0))

    ax2 = axes[0, 1]
    cc = b_df.cohort.value_counts().reindex(cohorts).fillna(0).sort_values()
    ax2.barh(range(len(cc)), cc.values, color=C_SIG, alpha=0.85, height=0.8, zorder=3)
    ax2.set_yticks(range(len(cc))); ax2.set_yticklabels(cc.index, fontsize=FS_TICK)
    axis_label(ax2, x="Number of pairs (IG gain ≥ 50%)")
    style_axes(ax2, grid="x")

    ax3 = axes[1, 0]
    cohorts_c = sorted(c_df.cohort.unique())
    lo_, hi_ = np.nanpercentile(c_df.log_gain.dropna(), [0.5, 99.5])
    xs = np.linspace(lo_, hi_, 300)
    for i, co in enumerate(cohorts_c):
        v = c_df.loc[c_df.cohort == co, "log_gain"].dropna().values
        if len(v) < 2:
            continue
        dens = gaussian_kde(v, bw_method=0.3)(xs)
        dens = dens / (dens.max() + 1e-9) * 0.8
        col = pal[cohorts.index(co)]
        ax3.fill_between(xs, i, i + dens, alpha=0.6, color=col, lw=0)
        ax3.plot(xs, i + dens, color=col, lw=0.6)
    ax3.set_xlim(lo_, hi_)
    ax3.set_yticks(np.arange(len(cohorts_c)) + 0.4); ax3.set_yticklabels(cohorts_c, fontsize=FS_TICK)
    axis_label(ax3, x="log$_{10}$ IG gain % (total IG ≥ 5)", y="Cohort")
    style_axes(ax3, grid=None)

    ax4 = axes[1, 1]
    sns.violinplot(data=s, x="cohort", y="pct_med", order=cohorts, hue="cohort", hue_order=cohorts, palette=pal,
                   legend=False, ax=ax4, width=0.8, linewidth=0.6, inner_kws={"box_width": 2, "whis_width": 0.5})
    ax4.axhline(100, color="gray", ls="--", lw=0.6)
    axis_label(ax4, x="Cohort", y="Total IG (% of cohort median)")
    plt.setp(ax4.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=FS_TICK)
    style_axes(ax4, grid="y")
    labels_row(fig, [(ax, "a"), (ax2, "b")])
    labels_row(fig, [(ax3, "c"), (ax4, "d")])
    save(fig, out, "Figure_4", formats, dpi)


# ------------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="high_res")
    ap.add_argument("--only", default="1,2,3,4", help="comma separated figure numbers")
    ap.add_argument("--formats", default="pdf,png")
    ap.add_argument("--dpi", type=int, default=600)
    ap.add_argument("--fig1-data")
    ap.add_argument("--fig2-data")
    ap.add_argument("--fig3-candidates")
    ap.add_argument("--prepared-data")
    ap.add_argument("--fig3-cache")
    ap.add_argument("--fig4-data")
    ap.add_argument("--no-arrows", action="store_true", help="Figure 4: omit the two arrows on panel a")
    a = ap.parse_args()
    fm = a.formats.split(",")
    todo = {int(i) for i in a.only.split(",")}
    print("font:", plt.rcParams["font.family"])
    if 1 in todo:
        print("Figure 1"); figure1(a.fig1_data, a.out, fm, a.dpi)
    if 2 in todo:
        print("Figure 2"); figure2(a.fig2_data, a.out, fm, a.dpi)
    if 3 in todo:
        print("Figure 3")
        figure3(a.fig3_candidates, a.prepared_data, a.fig3_cache or os.path.join(a.out, "figure3_effects.csv"), a.out, fm, a.dpi)
    if 4 in todo:
        print("Figure 4"); figure4(a.fig4_data, a.out, fm, a.dpi, arrows=not a.no_arrows)


if __name__ == "__main__":
    main()
