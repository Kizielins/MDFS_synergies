#!/usr/bin/env python3
"""Figure 2 sandbox (old Figure 3): hub-pair selection frequency + per-cohort effect sizes.

Same pair selection as R2_classification_performance/effect_sizes.ipynb (Figure 3 cell), but reads the saved
tables instead of recomputing the Mann-Whitney tests.

usage: python make_figure2.py <dir with cohort_level_tests.csv, final_synthetic_features_table.csv,
                              meta_results_all_pairs.csv> [outdir]
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(__file__))

INK, MUTED, GRID = "#1B2330", "#6A707A", "#E8EBEF"
C_BASE, C_CONTRIB, C_PAIR = "#B9C7DD", "#F2C777", "#3E9B77"
C_SIG, C_NS = "#1B6A86", "#C5C9D0"
plt.rcParams.update({"font.family": ["Helvetica", "Arial", "DejaVu Sans"], "svg.fonttype": "none",
                     "mathtext.fontset": "dejavusans"})

tests = pd.read_csv(os.path.join(DATA, "cohort_level_tests.csv"))
sel = pd.read_csv(os.path.join(DATA, "final_synthetic_features_table.csv"))
meta = pd.read_csv(os.path.join(DATA, "meta_results_all_pairs.csv"))

ranked = meta.merge(sel[["Data Type", "Feature", "Selection % Gain"]].rename(
    columns={"Data Type": "data_type", "Feature": "feature_pair"}), on=["data_type", "feature_pair"], how="left")
pm = ranked[ranked.feature_pair.str.contains("Parvimonas_micra")]
top = pm[pm["Selection % Gain"] >= 99.0]
r1 = (top if len(top) else pm).loc[(top if len(top) else pm).z_meta.abs().idxmax()]
fn = ranked[ranked.feature_pair.str.contains("Fusobacterium_nucleatum") & (ranked.feature_pair != r1.feature_pair)]
r2 = fn.loc[fn.z_meta.abs().idxmax()]
pairs = [(r1.data_type, r1.feature_pair), (r2.data_type, r2.feature_pair)]


def nice(x):
    return x.replace("_", " ")


def cohort_label(c):
    c = c.replace("Public_study__", "").replace("This_study_", "")
    if c.startswith("cohort"):
        c = c.split("__", 1)[1].replace("_", " ")
    return c.replace("_", " ")


def italic_species(x):
    return "$\\it{" + nice(x).replace(" ", "\\ ") + "}$"


fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.4), gridspec_kw={"width_ratios": [1.2, 1.1], "hspace": 0.62,
                                                              "wspace": 0.5})
for row, (dt, fp) in enumerate(pairs):
    pd_ = tests[(tests.data_type == dt) & (tests.feature_pair == fp)].copy()
    s = sel[(sel["Data Type"] == dt) & (sel.Feature == fp)].iloc[0]
    f = [(s.Feature1, s["Feature1_Selection%"]), (s.Feature2, s["Feature2_Selection%"])]
    f.sort(key=lambda t: t[1])
    (bn, bv), (cn, cv) = f
    pv = s["Selection %"]
    avg = np.mean([bv, cv])

    ax = axes[row, 0]
    ax.bar([0, 1, 2], [bv, cv, pv], width=0.62, color=[C_BASE, C_CONTRIB, C_PAIR], zorder=3)
    ax.axhline(avg, color=INK, ls="--", lw=0.9, zorder=4, label=f"mean of single features ({avg:.0f}%)")
    ax.legend(loc="upper left", frameon=False, fontsize=8, labelcolor=MUTED, handlelength=1.6)
    for x, v in enumerate([bv, cv, pv]):
        ax.text(x, v + 1.5, f"{v:.0f}%", ha="center", va="bottom", fontsize=9, color=INK, fontweight="bold")
    ax.annotate("", xy=(2.43, pv), xytext=(2.43, avg), arrowprops=dict(arrowstyle="<->", lw=0.9, color=INK))
    ax.text(2.5, (pv + avg) / 2, f"+{pv - avg:.0f} pp", ha="left", va="center", fontsize=9, fontweight="bold",
            color=C_PAIR)
    ax.set_xlim(-0.5, 3.1); ax.set_ylim(0, 125); ax.set_yticks(range(0, 101, 20))
    ax.set_xticks([0, 1, 2])
    def tick(n, role):
        g, sp = nice(n).split(" ", 1)
        return f"$\\it{{{g}}}$\n$\\it{{{sp}}}$\n({role})"
    ax.set_xticklabels([tick(bn, "base"), tick(cn, "contributor"), "Synthetic\npair"], fontsize=7.5)
    ax.set_ylabel("Selection frequency (% of LOCO folds)", fontsize=9, color=MUTED)
    ax.set_title("Selection frequency", loc="left", fontsize=10.5, fontweight="bold", color=INK, pad=8)

    ax2 = axes[row, 1]
    pd_ = pd_.assign(name=pd_.cohort.map(cohort_label)).sort_values("rbc", ascending=True)
    y = np.arange(len(pd_))
    sig = pd_.p_value_fdr.fillna(1) < 0.05
    ax2.barh(y, pd_.rbc, height=0.62, color=np.where(sig, C_SIG, C_NS), zorder=3)
    ax2.axvline(0, color=INK, lw=0.9, zorder=4)
    ax2.set_yticks(y); ax2.set_yticklabels(pd_.name, fontsize=8)
    ax2.set_xlabel("Rank-biserial r (CRC vs control)", fontsize=9, color=MUTED)
    n_sig = int(sig.sum()); n_cons = int((pd_.rbc > 0).sum() if pd_.rbc.median() > 0 else (pd_.rbc < 0).sum())
    ax2.set_title("Effect size per cohort", loc="left", fontsize=10.5, fontweight="bold", color=INK, pad=22)
    ax2.text(0, 1.03, f"{n_sig}/{len(pd_)} cohorts significant (FDR q < 0.05), {n_cons}/{len(pd_)} same direction",
             transform=ax2.transAxes, fontsize=8, color=MUTED)
    for a in (ax, ax2):
        a.grid(axis="x" if a is ax2 else "y", color=GRID, lw=0.6, zorder=0); a.set_axisbelow(True)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            a.spines[sp].set_color("#9AA1AB")
        a.tick_params(colors=MUTED, labelsize=8, length=3)
    ax2.legend(handles=[Patch(fc=C_SIG, label="FDR q < 0.05"),
                        Patch(fc=C_NS, label="FDR q \u2265 0.05")],
               loc="lower right", frameon=False, fontsize=8, labelcolor=MUTED, handlelength=1.2)
    if row == 0:
        ax.text(-0.30, 1.14, "a", transform=ax.transAxes, fontsize=16, fontweight="bold", color=INK, va="top")
        ax2.text(-0.62, 1.14, "b", transform=ax2.transAxes, fontsize=16, fontweight="bold", color=INK, va="top")
    def sp(n):
        g, c = nice(n).split(" ", 1)
        return f"$\\it{{{g[0]}.\\ {c}}}$"
    ax.text(0.0, 1.30, f"Pair {row + 1}:  {sp(bn)} (base) + {sp(cn)} (contributor)", transform=ax.transAxes,
            fontsize=10, fontweight="bold", color=INK, va="top")
    print(row, dt, fp, f"base={bv} contrib={cv} pair={pv} sig={n_sig}/{len(pd_)}")

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"figure2.{ext}"), dpi=200, bbox_inches="tight", facecolor="white")
print("saved")
