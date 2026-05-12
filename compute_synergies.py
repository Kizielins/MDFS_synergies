#!/usr/bin/env python3
"""
Compute metagenomic feature pairings and synergies via MDFS information gain.

For every pair of features identified by MDFS 2D, reports:
  - base_IG_1D: the base feature's individual 1D information gain
  - contributing_IG_1D: the contributing feature's individual 1D information gain
  - IG_2D_added: additional IG from the contributing feature conditioned on the base
  - total_IG: joint IG of the pair (base_IG_1D + IG_2D_added)
  - ig_gain_pct: IG_2D_added / base_IG_1D * 100

MDFS is run N_RUNS times with different seeds and IGs are averaged across runs for
robustness.  The output table is sorted by ig_gain_pct (descending).

Usage:
    python compute_synergies.py --X samples_features.tsv --y labels.tsv
    python compute_synergies.py --X samples_features.tsv --y labels.tsv --out results.txt
    python compute_synergies.py --X samples_features.tsv --y labels.tsv --n-runs 5

Input formats:
    X: tab-separated, first column = sample IDs, remaining columns = feature values
    y: tab-separated, first column = sample IDs, second column = labels (0=disease, 1=control)

Output (tab-separated, stdout or --out file):
    #  base  contributing  base_IG_1D  contributing_IG_1D  IG_2D_added  total_IG  ig_gain_pct
"""

import argparse
import os
import subprocess
import sys
import tempfile

import pandas as pd

N_RUNS_DEFAULT = 3
_SEEDS = [12, 13, 14, 15, 16]

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def _find_run_mdfs_r():
    candidates = [
        os.path.join(SCRIPT_DIR, "run_mdfs.R"),
        os.path.join(SCRIPT_DIR, "run_mdfs.r"),
        "run_mdfs.R",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return os.path.abspath(c)
    raise FileNotFoundError(
        "run_mdfs.R not found next to this script. "
        "Ensure run_mdfs.R is in the same directory as compute_synergies.py."
    )


def _run_mdfs_r(X_path, y_path, out_2d, out_unique, out_1d, seed):
    rscript = _find_run_mdfs_r()
    cmd = ["Rscript", rscript, X_path, y_path, out_2d, out_unique, out_1d, str(seed)]
    result = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"run_mdfs.R failed (exit {result.returncode}).\n"
            f"STDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )


def _one_mdfs_run(X: pd.DataFrame, y: pd.Series, tmp_dir: str, seed: int):
    """Run MDFS once; return pairs_df with columns: base, contributing,
    base_IG_1D, contributing_IG_1D, IG_2D_added, total_IG."""
    X_path    = os.path.join(tmp_dir, f"X_{seed}.tsv")
    y_path    = os.path.join(tmp_dir, f"y_{seed}.tsv")
    out_2d    = os.path.join(tmp_dir, f"2d_{seed}.tsv")
    out_uniq  = os.path.join(tmp_dir, f"uniq_{seed}.txt")
    out_1d    = os.path.join(tmp_dir, f"1d_{seed}.tsv")

    X.to_csv(X_path, sep="\t")
    y.rename("label").to_csv(y_path, sep="\t", header=True)

    _run_mdfs_r(X_path, y_path, out_2d, out_uniq, out_1d, seed)

    # --- 2D pairs (new format from run_mdfs.R) ---
    pairs = pd.read_csv(out_2d, sep="\t")
    if pairs.empty:
        return pd.DataFrame(columns=["base", "contributing", "base_IG_1D",
                                      "contributing_IG_1D", "IG_2D_added", "total_IG"])

    pairs["base"] = pairs["base"].astype(str).str.replace(" ", "_")
    pairs["contributing"] = pairs["contributing"].astype(str).str.replace(" ", "_")
    for col in ["base_IG_1D", "contributing_IG_1D", "IG_2D_added", "total_IG"]:
        pairs[col] = pairs[col].astype(float)

    return pairs[["base", "contributing", "base_IG_1D", "contributing_IG_1D",
                   "IG_2D_added", "total_IG"]]


def compute_synergies(X: pd.DataFrame, y: pd.Series, n_runs: int = N_RUNS_DEFAULT) -> pd.DataFrame:
    """
    Run MDFS n_runs times (averaged) and compute synergy metrics for all feature pairs.

    Returns a DataFrame with columns:
        base, contributing, base_IG_1D, contributing_IG_1D, IG_2D_added, total_IG, ig_gain_pct
    sorted by ig_gain_pct descending.
    """
    seeds = _SEEDS[:n_runs]
    all_rows = []

    with tempfile.TemporaryDirectory(prefix="mdfs_syn_") as tmp:
        for i, seed in enumerate(seeds):
            print(f"  Run {i + 1}/{n_runs} (seed={seed})...", file=sys.stderr)
            pairs = _one_mdfs_run(X, y, tmp, seed)
            if not pairs.empty:
                pairs["_seed"] = seed
                all_rows.append(pairs)

    if not all_rows:
        return pd.DataFrame(columns=["base", "contributing", "base_IG_1D",
                                      "contributing_IG_1D", "IG_2D_added",
                                      "total_IG", "ig_gain_pct"])

    info = pd.concat(all_rows, ignore_index=True)

    # Average across runs
    stats = info.groupby(["base", "contributing"]).agg(
        base_IG_1D=("base_IG_1D", "mean"),
        contributing_IG_1D=("contributing_IG_1D", "mean"),
        IG_2D_added=("IG_2D_added", "mean"),
        total_IG=("total_IG", "mean"),
    ).reset_index()

    # ig_gain_pct = IG_2D_added / base_IG_1D * 100
    denom = stats["base_IG_1D"].replace(0.0, float("nan"))
    stats["ig_gain_pct"] = (stats["IG_2D_added"] / denom) * 100

    stats = stats.sort_values("ig_gain_pct", ascending=False).reset_index(drop=True)
    return stats


def write_table(stats: pd.DataFrame, out_path: str = None) -> None:
    """Write tab-separated table with row numbers to out_path or stdout."""
    header = "#\tbase\tcontributing\tbase_IG_1D\tcontributing_IG_1D\tIG_2D_added\ttotal_IG\tig_gain_pct"
    lines = [header]
    for i, row in stats.iterrows():
        lines.append(
            f"{i + 1}\t{row['base']}\t{row['contributing']}\t"
            f"{row['base_IG_1D']:.6f}\t{row['contributing_IG_1D']:.6f}\t"
            f"{row['IG_2D_added']:.6f}\t{row['total_IG']:.6f}\t"
            f"{row['ig_gain_pct']:.2f}"
        )
    content = "\n".join(lines) + "\n"
    if out_path:
        with open(out_path, "w") as fh:
            fh.write(content)
        print(f"Wrote {len(stats)} pairs to {out_path}", file=sys.stderr)
    else:
        sys.stdout.write(content)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute metagenomic feature pairings / synergies using MDFS information gain.\n"
            "Outputs a tab-separated table sorted by ig_gain_pct (descending)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--X", required=True,
        help="Feature matrix TSV (samples x features; first column = sample IDs)",
    )
    parser.add_argument(
        "--y", required=True,
        help="Labels TSV (two columns: sample_id, label; 0 = disease, 1 = control/healthy)",
    )
    parser.add_argument(
        "--out", default=None,
        help="Output file path (default: stdout)",
    )
    parser.add_argument(
        "--n-runs", type=int, default=N_RUNS_DEFAULT,
        help=f"Number of MDFS runs to average over (default: {N_RUNS_DEFAULT})",
    )
    args = parser.parse_args()

    X = pd.read_csv(args.X, sep="\t", index_col=0)
    y = pd.read_csv(args.y, sep="\t", index_col=0).squeeze()
    y = y.astype(int)

    common = X.index.intersection(y.index)
    if len(common) == 0:
        print("Error: no shared sample IDs between X and y.", file=sys.stderr)
        sys.exit(1)
    X = X.loc[common]
    y = y.loc[common]

    print(f"Loaded {len(X)} samples, {X.shape[1]} features.", file=sys.stderr)

    stats = compute_synergies(X, y, n_runs=args.n_runs)
    if stats.empty:
        print("No pairs found.", file=sys.stderr)
        sys.exit(1)

    write_table(stats, args.out)


if __name__ == "__main__":
    main()
