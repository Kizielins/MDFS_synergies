#!/usr/bin/env python3
"""
Compute bacterial feature pairings and synergies via MDFS information gain.

For every pair of features identified by MDFS, reports:
  - 1D IG for each feature individually
  - 2D IG for the pair jointly
  - % synergy gain = (pair_ig - median(f1_ig, f2_ig)) / median(f1_ig, f2_ig) * 100

MDFS is run N_RUNS times with different seeds and IGs are averaged across runs for
robustness.  The output table is sorted by % synergy gain (descending).

Usage:
    python compute_synergies.py --X samples_features.tsv --y labels.tsv
    python compute_synergies.py --X samples_features.tsv --y labels.tsv --out results.txt
    python compute_synergies.py --X samples_features.tsv --y labels.tsv --n-runs 5

Input formats:
    X: tab-separated, first column = sample IDs, remaining columns = feature values
    y: tab-separated, first column = sample IDs, second column = labels (0=disease, 1=control)

Output (tab-separated, stdout or --out file):
    #  Feature_1  Feature_2  IG_F1  IG_F2  IG_Pair  Pct_Synergy_Gain
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
    """Run MDFS once; return (pairs_df, ig_1d_dict)."""
    X_path    = os.path.join(tmp_dir, f"X_{seed}.tsv")
    y_path    = os.path.join(tmp_dir, f"y_{seed}.tsv")
    out_2d    = os.path.join(tmp_dir, f"2d_{seed}.tsv")
    out_uniq  = os.path.join(tmp_dir, f"uniq_{seed}.txt")
    out_1d    = os.path.join(tmp_dir, f"1d_{seed}.tsv")

    X.to_csv(X_path, sep="\t")
    y.rename("label").to_csv(y_path, sep="\t", header=True)

    _run_mdfs_r(X_path, y_path, out_2d, out_uniq, out_1d, seed)

    # --- 2D pairs ---
    pairs = pd.read_csv(out_2d, sep="\t")
    if "Feature1" not in pairs.columns and pairs.shape[1] >= 3:
        pairs.columns = ["Feature1", "Feature2", "IG"] + list(pairs.columns[3:])
    pairs["Feature1"] = pairs["Feature1"].astype(str).str.replace(" ", "_")
    pairs["Feature2"] = pairs["Feature2"].astype(str).str.replace(" ", "_")

    ig_col = "IG" if "IG" in pairs.columns else pairs.columns[2]
    pairs[ig_col] = pairs[ig_col].astype(float)
    # Keep best IG per pair (MDFS may return multiple discretisations)
    pairs = pairs.groupby(["Feature1", "Feature2"], as_index=False)[ig_col].max()
    if ig_col != "IG":
        pairs = pairs.rename(columns={ig_col: "IG"})

    # --- 1D IGs ---
    ig_1d = pd.read_csv(out_1d, sep="\t")
    if "Feature" not in ig_1d.columns and ig_1d.shape[1] >= 2:
        ig_1d.columns = ["Feature", "IG"] + list(ig_1d.columns[2:])
    ig_1d["Feature"] = ig_1d["Feature"].astype(str).str.replace(" ", "_")
    ig_1d_dict = dict(zip(ig_1d["Feature"], ig_1d["IG"].astype(float)))

    return pairs, ig_1d_dict


def compute_synergies(X: pd.DataFrame, y: pd.Series, n_runs: int = N_RUNS_DEFAULT) -> pd.DataFrame:
    """
    Run MDFS n_runs times (averaged) and compute synergy metrics for all feature pairs.

    Returns a DataFrame with columns:
        feature1, feature2, f1_ig, f2_ig, pair_ig, pct_synergy_gain
    sorted by pct_synergy_gain descending.
    """
    seeds = _SEEDS[:n_runs]
    all_rows = []

    with tempfile.TemporaryDirectory(prefix="mdfs_syn_") as tmp:
        for i, seed in enumerate(seeds):
            print(f"  Run {i + 1}/{n_runs} (seed={seed})...", file=sys.stderr)
            pairs, ig_1d_dict = _one_mdfs_run(X, y, tmp, seed)
            for _, row in pairs.iterrows():
                f1, f2 = row["Feature1"], row["Feature2"]
                all_rows.append({
                    "feature1": f1,
                    "feature2": f2,
                    "pair_ig": float(row["IG"]),
                    "f1_ig": ig_1d_dict.get(f1, float("nan")),
                    "f2_ig": ig_1d_dict.get(f2, float("nan")),
                })

    if not all_rows:
        return pd.DataFrame(columns=["feature1", "feature2", "f1_ig", "f2_ig", "pair_ig", "pct_synergy_gain"])

    info = pd.DataFrame(all_rows)

    # Average across runs
    stats = info.groupby(["feature1", "feature2"]).agg(
        pair_ig=("pair_ig", "mean"),
        f1_ig=("f1_ig", "mean"),
        f2_ig=("f2_ig", "mean"),
    ).reset_index()

    # Median of the two individual IGs (= mean for 2 values)
    stats["median_individual_ig"] = stats[["f1_ig", "f2_ig"]].median(axis=1)

    # % synergy gain relative to median individual IG
    denom = stats["median_individual_ig"].replace(0.0, float("nan"))
    stats["pct_synergy_gain"] = (stats["pair_ig"] - stats["median_individual_ig"]) / denom * 100

    stats = stats.sort_values("pct_synergy_gain", ascending=False).reset_index(drop=True)
    return stats[["feature1", "feature2", "f1_ig", "f2_ig", "pair_ig", "pct_synergy_gain"]]


def write_table(stats: pd.DataFrame, out_path: str = None) -> None:
    """Write tab-separated table with row numbers to out_path or stdout."""
    header = "#\tFeature_1\tFeature_2\tIG_F1\tIG_F2\tIG_Pair\tPct_Synergy_Gain"
    lines = [header]
    for i, row in stats.iterrows():
        lines.append(
            f"{i + 1}\t{row['feature1']}\t{row['feature2']}\t"
            f"{row['f1_ig']:.6f}\t{row['f2_ig']:.6f}\t"
            f"{row['pair_ig']:.6f}\t{row['pct_synergy_gain']:.2f}"
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
            "Compute bacterial feature pairings / synergies using MDFS information gain.\n"
            "Outputs a tab-separated table sorted by % synergy gain (descending)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--X", required=True,
        help="Feature matrix TSV (samples × features; first column = sample IDs)",
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
