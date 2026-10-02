#!/usr/bin/env python3
"""Sensitivity of the pair search to the MDFS discretisation (Reviewer 1, comment 10).

LOCO fold-0 training set (taxa). Configurations, all with the pipeline's filters:
  default        divisions = 1 (two bins per feature), cut range = GetRange(n, dimensions = 2), seed 42
  default_seed43 same, seed 43 (reference for seed-to-seed variation)
  range_1D       cut range = GetRange(n, dimensions = 1) (wider range of cut points)
  range_half     cut range = half of the default
  divisions_2    divisions = 2 (three bins per feature; chi-square df = 2 for 1D and 6 for 2D)
For each configuration: number of 2D-significant features, pairs passing, and the overlap of the
top-2 % pair set with the default (Jaccard), plus the Spearman correlation of total IG of pairs
present in both. Output: discretisation_sensitivity.tsv
"""
import os
import re
import sys
import subprocess
import tempfile
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import training_set, unordered, top_fraction, HERE  # noqa: E402

RSCRIPT = r'''
suppressMessages({library(data.table); library(MDFS)})
args <- commandArgs(trailingOnly = TRUE)
xp <- args[1]; yp <- args[2]; out <- args[3]; seed <- as.integer(args[4]); divisions <- as.integer(args[5]); range_mode <- args[6]
X <- fread(xp, sep = "\t"); X <- as.matrix(X[, -1, with = FALSE]); y <- fread(yp, sep = "\t")[[2]]
n <- ncol(X); fn <- colnames(X)
rng <- GetRange(n = nrow(X), dimensions = 2, divisions = divisions)
if (range_mode == "range_1D") rng <- GetRange(n = nrow(X), dimensions = 1, divisions = divisions)
if (range_mode == "range_half") rng <- rng / 2
r1 <- MDFS(X, y, dimensions = 1, divisions = divisions, seed = seed, range = rng); ig1 <- r1$statistic
r2 <- MDFS(X, y, dimensions = 2, divisions = divisions, seed = seed, range = rng); sig <- r2$relevant.variables
cat("2D significant features:", length(sig), "\n")
df2 <- 1 * divisions * (divisions + 1)^(2 - 1)
thr <- qchisq(1 - 0.05 / max(1, n * (n - 1) / 2), df = df2) / (2 * log(2))
cat("threshold", thr, "\n")
res <- data.frame()
if (length(sig) > 0) {
  tp <- ComputeInterestingTuples(X, y, dimensions = 2, divisions = divisions, seed = seed, range = rng,
                                 interesting.vars = sig, I.lower = ig1, ig.thr = thr)
  if (!is.null(tp) && nrow(tp) > 0) {
    base <- ifelse(tp$Var == tp$Tuple.1, tp$Tuple.2, tp$Tuple.1)
    res <- data.frame(base = fn[base], contributing = fn[tp$Var], base_IG_1D = ig1[base], contributing_IG_1D = ig1[tp$Var],
                      IG_2D_added = tp$IG, total_IG = ig1[base] + tp$IG)
    res <- res[res$total_IG > pmax(res$base_IG_1D, res$contributing_IG_1D), ]
  }
}
write.table(res, out, sep = "\t", row.names = FALSE)
'''
CONFIGS = [("default", 42, 1, "default"), ("default_seed43", 43, 1, "default"), ("range_1D", 42, 1, "range_1D"),
           ("range_half", 42, 1, "range_half"), ("divisions_2", 42, 2, "default")]


def main():
    X, y = training_set("taxa")
    rows, sets, tables = [], {}, {}
    with tempfile.TemporaryDirectory(prefix="r5disc_") as tmp:
        rp = os.path.join(tmp, "disc.R"); open(rp, "w").write(RSCRIPT)
        xp, yp = os.path.join(tmp, "X.tsv"), os.path.join(tmp, "y.tsv")
        X.to_csv(xp, sep="\t"); y.rename("label").to_csv(yp, sep="\t", header=True)
        for name, seed, div, mode in CONFIGS:
            op = os.path.join(tmp, f"{name}.tsv")
            r = subprocess.run(["Rscript", rp, xp, yp, op, str(seed), str(div), mode], capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(name + r.stderr[-800:])
            d = pd.read_csv(op, sep="\t")
            sets[name] = (unordered(d), top_fraction(d)); tables[name] = d
            rows.append({"config": name, "seed": seed, "divisions": div, "range": mode,
                         "features_2D_significant": int(re.search(r"2D significant features: (\d+)", r.stdout).group(1)),
                         "pair_IG_threshold_bits": float(re.search(r"threshold ([\d.]+)", r.stdout).group(1)),
                         "unordered_pairs_passing": len(sets[name][0]), "top2pct_pairs": len(sets[name][1])})
            print(rows[-1], flush=True)
    ref_all, ref_top = sets["default"]
    dref = tables["default"].copy(); dref["key"] = [tuple(sorted((a, b))) for a, b in zip(dref.base, dref.contributing)]
    dref = dref.groupby("key").total_IG.max()
    for r in rows:
        a, t = sets[r["config"]]
        r["jaccard_all_vs_default"] = round(len(a & ref_all) / len(a | ref_all), 3) if (a | ref_all) else None
        r["jaccard_top2pct_vs_default"] = round(len(t & ref_top) / len(t | ref_top), 3) if (t | ref_top) else None
        d = tables[r["config"]].copy(); d["key"] = [tuple(sorted((x, z))) for x, z in zip(d.base, d.contributing)]
        d = d.groupby("key").total_IG.max(); common = d.index.intersection(dref.index)
        r["spearman_total_IG_shared_pairs"] = round(spearmanr(d.loc[common], dref.loc[common]).correlation, 3) if len(common) > 2 else None
        r["n_shared_pairs"] = len(common)
    out = pd.DataFrame(rows); out.to_csv(os.path.join(HERE, "discretisation_sensitivity.tsv"), sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
