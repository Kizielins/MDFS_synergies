#!/usr/bin/env python3
"""Which features do different selectors find, and does MDFS-2D cover them? (Reviewer 1, comment 2; Reviewer 2, comment 1)

For every LOCO fold (pooled training cohorts) and feature space, the set of relevant features is determined with
  MDFS-1D, MDFS-2D   relevant in >= 6 of 10 MDFS runs (seeds 42-51; Holm 0.05; same discretisation as run_mdfs.R)
  t-test             Welch t-test on log10(x + 1e-6), Holm 0.05
  U-test             Mann-Whitney U test, Holm 0.05
  LASSO              L1 logistic regression on standardised log10(x + 1e-6); penalty from 9 values (1e-3..10) by 5-fold CV (AUC),
                     liblinear, max 300 iterations; non-zero coefficients
  ReliefF            k = 10 nearest hits / misses, Manhattan distance on min-max scaled log10(x + 1e-6); it has no
                     significance test, so the top-N features are taken with N = number of MDFS-2D features
Coverage[A, B] = share of the features selected by A that B also selects, averaged over the 14 folds.
Outputs: method_coverage_<dt>.tsv, method_coverage_sizes.tsv, mdfs2d_only_features_<dt>.tsv, r5_numbers_coverage.tsv,
         method_coverage.png and ../manuscript/SuppFigure3_method_coverage.{png,pdf}
"""
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind
from sklearn.linear_model import LogisticRegressionCV
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import training_set, HERE, REPO  # noqa: E402

SEEDS = list(range(42, 52)); MIN_RUNS = 6
METHODS = ["MDFS-2D", "MDFS-1D", "t-test", "U-test", "LASSO", "ReliefF"]
RSCRIPT = r'''
suppressMessages({library(data.table); library(MDFS)})
a <- commandArgs(trailingOnly = TRUE); X <- fread(a[1], sep = "\t"); X <- as.matrix(X[, -1, with = FALSE]); y <- fread(a[2], sep = "\t")[[2]]
rng <- GetRange(n = nrow(X), dimensions = 2); n <- ncol(X); r1 <- integer(n); r2 <- integer(n); ig <- numeric(n)
for (seed in as.integer(strsplit(a[4], ",")[[1]])) {
  m1 <- MDFS(X, y, dimensions = 1, seed = seed, range = rng); m2 <- MDFS(X, y, dimensions = 2, seed = seed, range = rng)
  r1[m1$relevant.variables] <- r1[m1$relevant.variables] + 1L; r2[m2$relevant.variables] <- r2[m2$relevant.variables] + 1L; ig <- ig + m1$statistic
}
fwrite(data.table(feature = colnames(X), n_rel_1d = r1, n_rel_2d = r2, ig_1d = ig / length(strsplit(a[4], ",")[[1]])), a[3], sep = "\t")
'''


def relieff(Z, y, k=10):
    Z = (Z - Z.min(0)) / np.where(Z.max(0) > Z.min(0), Z.max(0) - Z.min(0), 1.0)
    w = np.zeros(Z.shape[1])
    for c in (0, 1):
        own, other = Z[y == c], Z[y != c]
        hit = NearestNeighbors(n_neighbors=k + 1, metric="manhattan").fit(own).kneighbors(own, return_distance=False)[:, 1:]
        mis = NearestNeighbors(n_neighbors=k, metric="manhattan").fit(other).kneighbors(own, return_distance=False)
        for i in range(len(own)):
            w += np.abs(own[i] - other[mis[i]]).mean(0) - np.abs(own[i] - own[hit[i]]).mean(0)
    return w / len(Z)


CACHE = os.path.join(HERE, "method_coverage_cache")


def one_fold(dt, fold, test_cohort):
    cp = os.path.join(CACHE, f"{dt}_{fold}.json")
    if os.path.exists(cp):
        c = json.load(open(cp)); return dt, fold, {k: set(v) for k, v in c["sel"].items()}, c["ig"]
    X, y = training_set(dt, test_cohort)
    with tempfile.TemporaryDirectory(prefix="cov_") as tmp:
        rp, xp, yp, op = (os.path.join(tmp, n) for n in ("m.R", "X.tsv", "y.tsv", "o.tsv"))
        open(rp, "w").write(RSCRIPT); X.to_csv(xp, sep="\t"); y.rename("label").to_csv(yp, sep="\t", header=True)
        r = subprocess.run(["Rscript", rp, xp, yp, op, ",".join(map(str, SEEDS))], capture_output=True, text=True)
        if r.returncode != 0: raise RuntimeError(r.stderr[-600:])
        m = pd.read_csv(op, sep="\t")
    feats = np.array(X.columns); yv = y.values; L = np.log10(X.values + 1e-6)
    sel = {"MDFS-1D": set(m.feature[m.n_rel_1d >= MIN_RUNS]), "MDFS-2D": set(m.feature[m.n_rel_2d >= MIN_RUNS])}
    pt = ttest_ind(L[yv == 0], L[yv == 1], equal_var=False).pvalue; pt = np.where(np.isfinite(pt), pt, 1.0)
    sel["t-test"] = set(feats[multipletests(pt, alpha=0.05, method="holm")[0]])
    pu = np.array([mannwhitneyu(X.values[yv == 0, j], X.values[yv == 1, j], alternative="two-sided").pvalue if X.values[:, j].std() > 0 else 1.0 for j in range(X.shape[1])])
    sel["U-test"] = set(feats[multipletests(pu, alpha=0.05, method="holm")[0]])
    Zs = StandardScaler().fit_transform(L)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lasso = LogisticRegressionCV(Cs=np.logspace(-3, 1, 9), cv=5, penalty="l1", solver="liblinear", scoring="roc_auc", class_weight="balanced",
                                     random_state=42, max_iter=300, tol=1e-3).fit(Zs, yv)
    sel["LASSO"] = set(feats[np.abs(lasso.coef_[0]) > 0])
    w = relieff(L.copy(), yv); n2 = max(1, len(sel["MDFS-2D"]))
    sel["ReliefF"] = set(feats[np.argsort(-w)[:n2]])
    ig = dict(zip(m.feature, m.ig_1d)); os.makedirs(CACHE, exist_ok=True)
    json.dump({"sel": {k: sorted(v) for k, v in sel.items()}, "ig": ig}, open(cp, "w"))
    return dt, fold, sel, ig


def main():
    tasks = []
    for dt in ("taxa", "functions", "combined"):
        f = pd.read_csv(os.path.join(REPO, "M2_loco_crc_pipeline", "results", dt, f"loco_results_{dt}.csv"))
        cohorts = [c for c in f["Test Cohort"] if c not in ("Mean", "Weighted Mean")]
        tasks += [(dt, i, c) for i, c in enumerate(cohorts)]
    res = {}
    with ProcessPoolExecutor(7) as ex:
        futs = [ex.submit(one_fold, *t) for t in tasks]
        for i, fu in enumerate(as_completed(futs), 1):
            dt, fold, sel, ig = fu.result(); res[(dt, fold)] = (sel, ig)
            print(f"[{i}/{len(tasks)}] {dt} fold {fold}: " + ", ".join(f"{k} {len(v)}" for k, v in sel.items()), flush=True)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    sizes, numbers = [], []
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.2))
    for ax, dt in zip(axes, ("taxa", "functions", "combined")):
        folds = sorted(f for d, f in res if d == dt)
        C = np.zeros((len(METHODS), len(METHODS)))
        for f in folds:
            sel = res[(dt, f)][0]
            for a, A in enumerate(METHODS):
                for b, B in enumerate(METHODS):
                    C[a, b] += (len(sel[A] & sel[B]) / len(sel[A])) if sel[A] else np.nan
        C /= len(folds)
        pd.DataFrame(C, index=METHODS, columns=METHODS).round(3).to_csv(os.path.join(HERE, f"method_coverage_{dt}.tsv"), sep="\t")
        n = {M: np.mean([len(res[(dt, f)][0][M]) for f in folds]) for M in METHODS}
        others = [M for M in METHODS if M != "MDFS-2D"]; uni = ["MDFS-1D", "t-test", "U-test"]
        only_all = [len(res[(dt, f)][0]["MDFS-2D"] - set().union(*[res[(dt, f)][0][M] for M in others])) for f in folds]
        only_uni = [len(res[(dt, f)][0]["MDFS-2D"] - set().union(*[res[(dt, f)][0][M] for M in uni])) for f in folds]
        miss = {M: np.mean([len(res[(dt, f)][0][M] - res[(dt, f)][0]["MDFS-2D"]) for f in folds]) for M in others}
        sizes.append({"data_type": dt, **{f"n_{M}": round(n[M], 1) for M in METHODS}, "MDFS-2D_only_vs_univariate": np.mean(only_uni), "MDFS-2D_only_vs_all": np.mean(only_all),
                      **{f"missed_by_MDFS-2D_{M}": round(miss[M], 1) for M in others}})
        cnt, igs = {}, {}
        for f in folds:
            sel, ig = res[(dt, f)]
            for x in sel["MDFS-2D"] - set().union(*[sel[M] for M in uni]):
                cnt[x] = cnt.get(x, 0) + 1; igs.setdefault(x, []).append(ig[x])
        pd.DataFrame([{"feature": k, "folds_MDFS-2D_only": v, "mean_1D_IG": np.mean(igs[k])} for k, v in cnt.items()]).sort_values("folds_MDFS-2D_only", ascending=False).to_csv(
            os.path.join(HERE, f"mdfs2d_only_features_{dt}.tsv"), sep="\t", index=False, float_format="%.2f")
        b2 = METHODS.index("MDFS-2D")
        numbers.append({"section": "R5", "claim": f"{dt}: share of each method's features also found by MDFS-2D (" + " / ".join(others) + ")",
                        "value": " / ".join(f"{100 * C[METHODS.index(M), b2]:.0f} %" for M in others), "source": f"R5 method_coverage_{dt}.tsv"})
        numbers.append({"section": "R5", "claim": f"{dt}: share of MDFS-2D features also found by each other method (" + " / ".join(others) + ")",
                        "value": " / ".join(f"{100 * C[b2, METHODS.index(M)]:.0f} %" for M in others), "source": f"R5 method_coverage_{dt}.tsv"})
        numbers.append({"section": "R5", "claim": f"{dt}: mean number of features per fold (" + " / ".join(METHODS) + "); MDFS-2D features found by no univariate method; by no other method",
                        "value": " / ".join(f"{n[M]:.0f}" for M in METHODS) + f"; {np.mean(only_uni):.1f}; {np.mean(only_all):.1f}", "source": "R5 method_coverage_sizes.tsv"})
        im = ax.imshow(C, vmin=0, vmax=1, cmap="viridis")
        ax.set_xticks(range(len(METHODS))); ax.set_xticklabels([f"{M}\n(n={n[M]:.0f})" for M in METHODS], fontsize=9); ax.set_yticks(range(len(METHODS))); ax.set_yticklabels(METHODS, fontsize=10)
        for a in range(len(METHODS)):
            for b in range(len(METHODS)):
                ax.text(b, a, f"{100 * C[a, b]:.0f}", ha="center", va="center", color="white" if C[a, b] < 0.6 else "black", fontsize=10)
        ax.set_title(dt.upper(), fontsize=14, fontweight="bold"); ax.set_xlabel("... are also selected by", fontsize=11)
        if dt == "taxa": ax.set_ylabel("% of the features selected by ...", fontsize=11)
    fig.colorbar(im, ax=axes, fraction=0.015, pad=0.01).set_label("share of row method's features covered by column method")
    for p in (os.path.join(HERE, "method_coverage.png"), os.path.join(REPO, "manuscript", "SuppFigure3_method_coverage.png"), os.path.join(REPO, "manuscript", "SuppFigure3_method_coverage.pdf")):
        fig.savefig(p, dpi=200, bbox_inches="tight")
    pd.DataFrame(sizes).to_csv(os.path.join(HERE, "method_coverage_sizes.tsv"), sep="\t", index=False, float_format="%.1f")
    pd.DataFrame(numbers).to_csv(os.path.join(HERE, "r5_numbers_coverage.tsv"), sep="\t", index=False)
    print(pd.DataFrame(sizes).to_string(index=False)); print(pd.DataFrame(numbers)[["claim", "value"]].to_string(index=False))


if __name__ == "__main__":
    main()
