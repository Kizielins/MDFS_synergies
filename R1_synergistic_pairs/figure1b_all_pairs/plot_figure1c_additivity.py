#!/usr/bin/env python3
"""
Figure 1c: is a pair more informative than its two members together?
One point per unordered pair, averaged over the LOCO folds x MDFS runs in which both features were present
(input: directed/<dt>_fold<NN>.tsv.gz written by run_all_directed.py).
  x = sum of the two individual (1D) IGs      y = total IG of the pair
  diagonal y = x: additive (interaction information = 0); above: super-additive; below: the two features overlap
  dashed line: diagonal + chance threshold (median over folds of the 99th percentile of the interaction information
               under within-class shuffling, R5_reviewer_checks/interaction_null_thresholds.tsv)
Colours: grey = never synergistic by the pipeline's definition; blue = synergistic in >= 1 fold-run;
         red = synergistic in >= 50 % of fold-runs AND mean interaction information above the chance threshold.
Writes figure1c_additivity.{png,pdf} here and ../../manuscript/Figure1c.{png,pdf}; an alternative view with the
excess over additivity on the y-axis (figure1c_alt_excess.png); and figure1c_summary.tsv.
"""
import glob
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MS = os.path.join(HERE, "..", "..", "manuscript")
THR = pd.read_csv(os.path.join(HERE, "..", "..", "R5_reviewer_checks", "interaction_null_thresholds.tsv"), sep="\t")
SHORT = {"Parvimonas_micra": "P. micra", "Peptostreptococcus_stomatis": "P. stomatis", "Fusobacterium_nucleatum": "F. nucleatum",
         "Dialister_pneumosintes": "D. pneumosintes", "Gemella_morbillorum": "G. morbillorum", "Streptococcus_mitis": "S. mitis",
         "Gemella_sanguinis": "G. sanguinis"}
LABELS = [("Parvimonas_micra", "Peptostreptococcus_stomatis"), ("Fusobacterium_nucleatum", "Peptostreptococcus_stomatis"),
          ("Dialister_pneumosintes", "Parvimonas_micra"), ("Gemella_morbillorum", "Streptococcus_mitis")]


def load(dt):
    a = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(os.path.join(HERE, "directed", f"{dt}_fold*.tsv.gz")))])
    a = a.groupby(["base", "contributing"], as_index=False).sum(numeric_only=True)
    a["ig_a"] = a.base_1D / a.n; a["tot"] = a.total / a.n; a["psyn"] = 100 * a.n_syn / a.n
    b = a.rename(columns={"base": "contributing", "contributing": "base", "ig_a": "ig_b", "psyn": "psyn_rev"})[["base", "contributing", "ig_b", "psyn_rev"]]
    m = a.merge(b, on=["base", "contributing"]); m = m[m.base < m.contributing].copy()
    m["sum1d"] = m.ig_a + m.ig_b; m["II"] = m.tot - m.sum1d; m["psyn"] = m[["psyn", "psyn_rev"]].max(axis=1)
    thr = THR[THR.data_type == dt].null_q99.median()
    m["cat"] = np.where(m.psyn == 0, "none", np.where((m.psyn >= 50) & (m.II > thr), "super", "syn"))
    return m, thr


data = {dt: load(dt) for dt in ("taxa", "functions", "combined")}
rows = []
for alt in (False, True):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.4), sharex=True, sharey=True)
    lim = max(max(d.sum1d.max(), d.tot.max()) for d, _ in data.values()) * 1.03
    for ax, (dt, (d, thr)) in zip(axes, data.items()):
        yv = (lambda q: q.II) if alt else (lambda q: q.tot)
        for cat, col, lab, al, sz in (("none", "#B0B0B0", "never synergistic", 0.35, 3), ("syn", "#3B7DD8", "synergistic (total IG above either feature)", 0.45, 3),
                                      ("super", "#C0392B", "super-additive (above chance threshold)", 0.8, 5)):
            q = d[d.cat == cat]
            ax.scatter(q.sum1d, yv(q), s=sz, c=col, alpha=al, linewidths=0, rasterized=True, label=f"{lab} (n={len(q):,})")
        if alt:
            ax.axhline(0, color="k", ls="--", lw=1, label="additive (total = sum)"); ax.axhline(thr, color="k", ls=":", lw=1, label=f"chance threshold (+{thr:.1f} bits)")
            ax.set_yscale("symlog", linthresh=10); ax.set_xlim(0, lim)
        else:
            ax.plot([0, lim], [0, lim], "k--", lw=1, label="additive (total = sum)"); ax.plot([0, lim], [thr, lim + thr], "k:", lw=1, label=f"chance threshold (+{thr:.1f} bits)")
            ax.set_xlim(0, lim); ax.set_ylim(0, lim * 0.75)
        for a, b in LABELS:
            q = d[(d.base == min(a, b)) & (d.contributing == max(a, b))]
            if len(q):
                ax.annotate(f"{SHORT.get(a, a)} + {SHORT.get(b, b)}", (q.sum1d.iloc[0], yv(q).iloc[0]), xytext=(-14, 14 if "micra" not in (a + b).lower() or "Dialister" in a else -16), textcoords="offset points",
                            fontsize=8, ha="right", arrowprops=dict(arrowstyle="-", lw=0.6))
        ax.set_title(f"{dt.upper()}  ({len(d):,} feature pairs)", fontsize=14, fontweight="bold"); ax.grid(alpha=0.3)
        ax.set_xlabel("Sum of the two individual (1D) IGs", fontsize=12); ax.legend(fontsize=8, loc="lower right" if not alt else "lower left", markerscale=3)
        if not alt:
            vc = d.cat.value_counts(); syn = d[d.psyn >= 50]
            rows.append({"data_type": dt, "pairs": len(d), "never_synergistic": int(vc.get("none", 0)), "synergistic": int(vc.get("syn", 0)), "super_additive": int(vc.get("super", 0)),
                         "chance_threshold_bits": round(thr, 2), "pairs_synergistic_in_ge_50pct_fold_runs": len(syn),
                         "share_of_those_super_additive": round((syn.II > thr).mean(), 4), "share_of_those_below_minus_threshold": round((syn.II < -thr).mean(), 4)})
    axes[0].set_ylabel("Total IG minus sum of individual IGs (interaction information)" if alt else "Total IG of the pair", fontsize=12)
    fig.tight_layout()
    if alt:
        fig.savefig(os.path.join(HERE, "figure1c_alt_excess.png"), dpi=200, bbox_inches="tight")
    else:
        for ext in ("png", "pdf"):
            fig.savefig(os.path.join(HERE, f"figure1c_additivity.{ext}"), dpi=200, bbox_inches="tight")
            fig.savefig(os.path.join(MS, f"Figure1c.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)
out = pd.DataFrame(rows); out.to_csv(os.path.join(HERE, "figure1c_summary.tsv"), sep="\t", index=False); print(out.to_string(index=False))
