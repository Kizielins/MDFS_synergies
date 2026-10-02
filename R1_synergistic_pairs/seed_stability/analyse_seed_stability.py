#!/usr/bin/env python3
"""
Summarise run_seed_stability.py output.

Per run: top-1% pair set (rank <= max(1, int(0.01 * n_unique_pairs))); the 2% grid is derived the same way.
Consensus(N, k): pairs in >= k of N runs.
Reproducibility: Jaccard between consensus sets built from two disjoint seed blocks
(seeds 42.. and 52..), averaged over folds.
Writes consensus_grid.tsv and selection_frequency.tsv.
"""
import glob
import os
from collections import Counter

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..")
BLOCK_A = list(range(42, 52))
BLOCK_B = list(range(52, 62))
GRID = {2: [1, 2], 3: [2, 3], 5: [3, 4, 5], 10: [5, 6, 7, 8, 10]}


def top1(path):
    d = pd.read_csv(path, sep="\t")
    n = max(1, int(d["n_unique_pairs"].iloc[0] * 0.01))
    return set(zip(d["feature_1"][:n], d["feature_2"][:n]))


def consensus(sets, k):
    c = Counter(p for s in sets for p in s)
    return {p for p, n in c.items() if n >= k}


def jac(a, b):
    return len(a & b) / len(a | b) if (a | b) else np.nan


rows, freq_rows = [], []
for dt in ["taxa", "functions", "combined"]:
    # optional comparison with the submitted (2-run) results, if those files are present
    _old_path = os.path.join(REPO, "archive", "submitted_results", f"mdfs_ig_info_{dt}.csv")
    old = pd.read_csv(_old_path) if os.path.exists(_old_path) else pd.DataFrame(columns=["fold", "run", "base", "contributing"])
    old = old[old["run"] == 0]
    folds = sorted({int(f.split("_fold")[1][:2]) for f in glob.glob(os.path.join(HERE, "pairs", f"{dt}_fold*"))})
    for fold in folds:
        runs = {s: top1(os.path.join(HERE, "pairs", f"{dt}_fold{fold:02d}_seed{s}.tsv.gz")) for s in BLOCK_A + BLOCK_B}
        o = old[old["fold"] == fold]
        old_set = {tuple(sorted(p)) for p in zip(o["base"], o["contributing"])}
        # single-run agreement
        pair_j = [jac(runs[a], runs[b]) for a, b in zip(BLOCK_A, BLOCK_B)]
        # selection frequency over all 20 seeds
        cnt = Counter(p for s in runs.values() for p in s)
        for p, n in cnt.items():
            freq_rows.append({"data_type": dt, "fold": fold, "freq": n / 20})
        for N, ks in GRID.items():
            A_sets = [runs[s] for s in BLOCK_A[:N]]
            B_sets = [runs[s] for s in BLOCK_B[:N]]
            for k in ks:
                A, B = consensus(A_sets, k), consensus(B_sets, k)
                rows.append({"data_type": dt, "fold": fold, "N": N, "k": k,
                             "top1_per_run": len(runs[42]), "single_run_jaccard": np.mean(pair_j),
                             "size": len(A), "jaccard_AB": jac(A, B),
                             "frac_in_old": len(A & old_set) / len(A) if A else np.nan,
                             "frac_old_kept": len(A & old_set) / len(old_set) if old_set else np.nan})

grid = pd.DataFrame(rows)
grid.to_csv(os.path.join(HERE, "consensus_grid.tsv"), sep="\t", index=False, float_format="%.4f")
freq = pd.DataFrame(freq_rows)
freq.to_csv(os.path.join(HERE, "selection_frequency.tsv"), sep="\t", index=False)

pd.set_option("display.width", 200)
print("Single-run agreement (Jaccard of top-1% sets, two seeds), mean over folds:")
print(grid.drop_duplicates(["data_type", "fold"]).groupby("data_type", sort=False)[["top1_per_run", "single_run_jaccard"]].mean().round(3), "\n")
summ = grid.groupby(["data_type", "N", "k"], sort=False).agg(
    size=("size", "mean"), jaccard_AB=("jaccard_AB", "mean"), jaccard_AB_min=("jaccard_AB", "min"),
    frac_in_old=("frac_in_old", "mean"), frac_old_kept=("frac_old_kept", "mean")).round(3)
print(summ.to_string(), "\n")
bins = [0, 0.2, 0.4, 0.6, 0.8, 0.999, 1.0]
freq["bin"] = pd.cut(freq["freq"], bins, labels=["<=20%", "20-40%", "40-60%", "60-80%", "80-<100%", "100%"])
print("Selection frequency over 20 seeds (share of pairs ever in a top-1% set):")
print(pd.crosstab(freq["data_type"], freq["bin"], normalize="index").round(3))
