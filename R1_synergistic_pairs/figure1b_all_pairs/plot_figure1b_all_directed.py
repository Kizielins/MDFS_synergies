#!/usr/bin/env python3
"""
Figure 1b, original methodology with all points.
One point per directed pair (base feature, contributing feature), averaged over the LOCO folds x MDFS
runs in which both features were present (up to 14 x 10).
x = base feature 1D IG;  y = base 1D IG + 2D IG added by the contributing feature.
Colour = % of fold-runs in which the pair passed significance and the synergy filter (grey = never). Writes figure1b_all_directed.png/.pdf and _summary.tsv.
"""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CATS = [  # drawn in this order
    ("not_sig", "not significant (added IG below threshold)", "#B0B0B0"),
    ("removed", "significant, removed by synergy filter", "#C0392B"),
    ("some", "synergistic in some runs", "#E0A800"),
    ("kept", "synergistic in every run", "#2E8B57"),
]


def load(dt):
    parts = [pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(HERE, "directed", f"{dt}_fold*.tsv.gz")))]
    a = pd.concat(parts).groupby(["base", "contributing"], as_index=False).sum(numeric_only=True)
    a["x"] = a.base_1D / a.n
    a["y"] = a.total / a.n
    a["cat"] = "not_sig"
    a.loc[a.n_sig > 0, "cat"] = "removed"
    a.loc[a.n_syn > 0, "cat"] = "some"
    a.loc[a.n_syn == a.n, "cat"] = "kept"
    return a


data = {dt: load(dt) for dt in ("taxa", "functions", "combined")}
lim = max(max(d.x.max(), d.y.max()) for d in data.values()) * 1.03
fig, axes = plt.subplots(1, 3, figsize=(18, 6.4), sharex=True, sharey=True)
rows = []
for ax, (dt, d) in zip(axes, data.items()):
    d["pct_syn"] = 100 * d.n_syn / d.n
    vc = d.cat.value_counts()
    rows.append({"data_type": dt, "directed_pairs": len(d), **{f"n_{c}": int(vc.get(c, 0)) for c, _, _ in CATS},
                 "n_syn_ge_50pct_runs": int((d.pct_syn >= 50).sum()),
                 "n_on_or_below_diagonal": int((d.y <= d.x).sum()),
                 "min_gap_y_minus_x": round((d.y - d.x).min(), 4)})
    never = d[d.n_syn == 0].sort_values("n_sig")
    ax.scatter(never.x, never.y, s=4, c="#B0B0B0", alpha=0.5, linewidths=0, rasterized=True,
               label=f"never synergistic (n={len(never):,})")
    syn = d[d.n_syn > 0].sort_values("pct_syn")          # most consistent drawn on top
    sc = ax.scatter(syn.x, syn.y, s=4, c=syn.pct_syn, cmap="viridis_r", vmin=0, vmax=100, alpha=0.6,
                    linewidths=0, rasterized=True, label=f"synergistic in >=1 run (n={len(syn):,})")
    ax.plot([0, lim], [0, lim], "k--", lw=1, label="y = x (no 2D gain)")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.grid(alpha=0.3)
    ax.set_title(f"{dt.upper()}  ({len(d):,} feature pairs)", fontsize=14, fontweight="bold")
    ax.set_xlabel("Base feature 1D IG", fontsize=12)
    ax.legend(fontsize=8, loc="upper left", markerscale=3)
cb = fig.colorbar(sc, ax=axes, fraction=0.015, pad=0.01)
cb.set_label("% of fold-runs in which the pair is synergistic (kept)", fontsize=11)
axes[0].set_ylabel("Base 1D IG + 2D IG added by the contributing feature", fontsize=12)
MS = os.path.join(HERE, "..", "..", "manuscript")
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(HERE, f"figure1b_all_directed.{ext}"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(MS, f"Figure1b.{ext}"), dpi=200, bbox_inches="tight")   # manuscript Figure 1b
pd.DataFrame(rows).to_csv(os.path.join(HERE, "figure1b_all_directed_summary.tsv"), sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))
