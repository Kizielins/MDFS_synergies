"""
MDFS synergistic pairs: run run_mdfs.R repeatedly and aggregate the runs.

One MDFS run (run_mdfs.R) returns every (base, contributing) pair that passes
the Bonferroni-corrected pair-level IG threshold and the true-synergy filter
total_IG > max(base_IG_1D, contributing_IG_1D). Each pair appears in both
directions.

Across runs (different seeds), a pair is counted as selected in a run when it
is among the top `top_fraction` of that run's unordered pairs, scored by the
larger total_IG of its two directions (at least one pair is kept per run).
Consensus pairs are those selected in at least `min_runs` runs.
"""

import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RUN_MDFS_R = os.path.join(SCRIPT_DIR, "run_mdfs.R")

PAIR_COLUMNS = ["base", "contributing", "base_IG_1D", "contributing_IG_1D", "IG_2D_added", "total_IG"]


def run_mdfs_once(X: pd.DataFrame, y: pd.Series, seed: int, tmp_dir: str) -> pd.DataFrame:
    """Run run_mdfs.R once and return its pairs (columns PAIR_COLUMNS; may be empty)."""
    paths = {k: os.path.join(tmp_dir, f"{k}_{seed}.tsv") for k in ("X", "y", "2d", "unique", "1d")}
    X.to_csv(paths["X"], sep="\t")
    y.rename("label").to_csv(paths["y"], sep="\t", header=True)

    cmd = ["Rscript", RUN_MDFS_R, paths["X"], paths["y"], paths["2d"], paths["unique"], paths["1d"], str(seed)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"run_mdfs.R failed (exit {result.returncode}).\n"
                           f"STDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}")

    try:
        pairs = pd.read_csv(paths["2d"], sep="\t")
    except pd.errors.EmptyDataError:
        pairs = pd.DataFrame()
    if pairs.empty or "base" not in pairs.columns:
        return pd.DataFrame(columns=PAIR_COLUMNS)
    return pairs[PAIR_COLUMNS]


def run_mdfs_repeated(X: pd.DataFrame, y: pd.Series, n_runs: int, seed: int, label: str = "") -> pd.DataFrame:
    """Run MDFS n_runs times with seeds seed, seed+1, ...; return all pairs with a `run` column (1-based)."""
    runs = []
    with tempfile.TemporaryDirectory(prefix="mdfs_") as tmp:
        for i in range(n_runs):
            print(f"    {label}MDFS run {i + 1}/{n_runs} (seed={seed + i})", file=sys.stderr)
            pairs = run_mdfs_once(X, y, seed + i, tmp)
            runs.append(pairs.assign(run=i + 1))
    return pd.concat(runs, ignore_index=True)


def _unordered(df: pd.DataFrame) -> pd.Series:
    """Unordered pair key (feature_1, feature_2) with feature_1 < feature_2."""
    return pd.Series([tuple(sorted(p)) for p in zip(df["base"], df["contributing"])], index=df.index, dtype=object)


def top_pairs_per_run(runs: pd.DataFrame, top_fraction: float) -> pd.DataFrame:
    """
    For each run, score unordered pairs by the larger total_IG of their two directions
    and keep the top `top_fraction` (at least one). Returns columns run, pair, score.
    """
    if runs.empty:
        return pd.DataFrame(columns=["run", "pair", "score"])
    scored = (runs.assign(pair=_unordered(runs))
                  .groupby(["run", "pair"])["total_IG"].max()
                  .rename("score").reset_index())
    kept = []
    for _, run_pairs in scored.groupby("run"):
        n_keep = max(1, int(len(run_pairs) * top_fraction))
        kept.append(run_pairs.nlargest(n_keep, "score"))
    return pd.concat(kept, ignore_index=True)


def consensus_pairs(runs: pd.DataFrame, top_fraction: float, min_runs: int) -> pd.DataFrame:
    """
    Pairs among the top `top_fraction` in at least `min_runs` runs.
    Columns: feature_1, feature_2, n_runs_top, total_IG_mean; sorted by n_runs_top, total_IG_mean.
    """
    cols = ["feature_1", "feature_2", "n_runs_top", "total_IG_mean"]
    top = top_pairs_per_run(runs, top_fraction)
    if top.empty:
        return pd.DataFrame(columns=cols)
    agg = top.groupby("pair").agg(n_runs_top=("run", "nunique"), total_IG_mean=("score", "mean")).reset_index()
    agg = agg[agg["n_runs_top"] >= min_runs]
    agg["feature_1"] = agg["pair"].str[0]
    agg["feature_2"] = agg["pair"].str[1]
    agg = agg.sort_values(["n_runs_top", "total_IG_mean"], ascending=False).reset_index(drop=True)
    return agg[cols]


def synergy_table(runs: pd.DataFrame, top_fraction: float, min_runs: int) -> pd.DataFrame:
    """
    One row per (base, contributing) direction, averaged over the runs in which the
    pair passed the MDFS filters, with:
      n_runs_significant  runs in which the pair passed both filters
      n_runs_top          runs in which the (unordered) pair was in the top fraction
      consensus           n_runs_top >= min_runs
    Sorted by total_IG_mean (descending).
    """
    cols = ["base", "contributing", "base_IG_1D_mean", "contributing_IG_1D_mean", "IG_2D_added_mean",
            "IG_2D_added_sd", "total_IG_mean", "ig_gain_pct", "n_runs_significant", "n_runs_top", "consensus"]
    if runs.empty:
        return pd.DataFrame(columns=cols)

    stats = runs.groupby(["base", "contributing"]).agg(
        base_IG_1D_mean=("base_IG_1D", "mean"),
        contributing_IG_1D_mean=("contributing_IG_1D", "mean"),
        IG_2D_added_mean=("IG_2D_added", "mean"),
        IG_2D_added_sd=("IG_2D_added", "std"),
        total_IG_mean=("total_IG", "mean"),
        n_runs_significant=("run", "nunique"),
    ).reset_index()
    stats["IG_2D_added_sd"] = stats["IG_2D_added_sd"].fillna(0.0)
    stats["ig_gain_pct"] = stats["IG_2D_added_mean"] / stats["base_IG_1D_mean"].replace(0.0, np.nan) * 100

    top = top_pairs_per_run(runs, top_fraction)
    n_top = top.groupby("pair")["run"].nunique() if not top.empty else pd.Series(dtype=int)
    stats["n_runs_top"] = _unordered(stats).map(n_top).fillna(0).astype(int)
    stats["consensus"] = stats["n_runs_top"] >= min_runs

    return stats.sort_values(["total_IG_mean", "base"], ascending=[False, True]).reset_index(drop=True)[cols]
