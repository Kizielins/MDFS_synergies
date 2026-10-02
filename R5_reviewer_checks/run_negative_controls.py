#!/usr/bin/env python3
"""Negative controls for the MDFS pair search (Reviewer 1, comment 7).

On the LOCO fold-0 training set of each data type, run_mdfs.R (Bonferroni pair threshold +
synergy filter, seed 42) is applied to
  real          the real labels (reference),
  label_perm    labels permuted across samples (global null: no feature carries disease signal),
  within_class  every feature shuffled independently within each class (marginal 1D signal of each
                feature is preserved, any dependence between features and hence any pairwise
                synergy is destroyed: null for the pair search specifically),
N_PERM times each. Output: negative_controls.tsv (one row per run) with the number of directed
pairs and unordered pairs passing the filters and the number of 2D-significant features.
"""
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import training_set, run_mdfs, unordered, HERE  # noqa: E402

N_PERM = 5
SEED = 42


def one(dt, kind, perm):
    X, y = training_set(dt)
    rng = np.random.default_rng(perm)
    if kind == "label_perm":
        y = pd.Series(rng.permutation(y.values), index=y.index)
    elif kind == "within_class":
        X = X.copy()
        for cls in (0, 1):
            idx = np.where(y.values == cls)[0]
            for j in range(X.shape[1]):
                X.iloc[idx, j] = X.iloc[rng.permutation(idx), j].values
    pairs, ig1d, out, wall, rss = run_mdfs(X, y, SEED)
    sig2d = int(re.search(r"2D significant features: (\d+)", out).group(1))
    sig1d = int(re.search(r"1D significant features: (\d+)", out).group(1))
    thr = re.search(r"Pair-level IG cutoff: ([\d.]+)", out)
    return {"data_type": dt, "control": kind, "permutation": perm, "n_samples": len(y), "n_features": X.shape[1],
            "features_1D_significant": sig1d, "features_2D_significant": sig2d,
            "pair_IG_threshold_bits": float(thr.group(1)) if thr else np.nan,
            "directed_pairs_passing": len(pairs), "unordered_pairs_passing": len(unordered(pairs)),
            "max_total_IG": float(pairs["total_IG"].max()) if len(pairs) else 0.0, "wall_s": round(wall)}


def main():
    tasks = [(dt, "real", 0) for dt in ("taxa", "functions", "combined")]
    tasks += [(dt, kind, p) for dt in ("taxa", "functions", "combined") for kind in ("label_perm", "within_class") for p in range(1, N_PERM + 1)]
    rows = []
    with ProcessPoolExecutor(4) as ex:
        futs = {ex.submit(one, *t): t for t in tasks}
        for i, fu in enumerate(as_completed(futs), 1):
            r = fu.result(); rows.append(r)
            print(f"[{i}/{len(tasks)}] {r['data_type']} {r['control']} {r['permutation']}: {r['unordered_pairs_passing']} pairs, "
                  f"{r['features_2D_significant']} 2D-significant features", flush=True)
    out = pd.DataFrame(rows).sort_values(["data_type", "control", "permutation"])
    out.to_csv(os.path.join(HERE, "negative_controls.tsv"), sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
