#!/usr/bin/env python3
"""
Run MDFS in 1D and 2D on cohort data and report information gains.

Runs MDFS 3 times per cohort (like the main pipeline), then aggregates IGs
across runs using mean and writes:
  - mdfs_2d_<cohort>.tsv: 2D tuples (base, contributing, base_IG_1D,
    contributing_IG_1D, IG_2D_added, total_IG) from last run
  - mdfs_1d_ig_<cohort>.tsv: per-feature 1D IG with p-values (from last run)
  - mdfs_ig_info_<cohort>.tsv: one row per run per (base, contributing) pair
  - ig_gain_summary_<cohort>.tsv: aggregated by (base, contributing) across runs

Usage:
  python run_mdfs_1d_2d.py --data-dir data --cohort VogtmannE_2016 [--out-dir results]
  python run_mdfs_1d_2d.py --data-dir data --all-cohorts [--out-dir results]
  python run_mdfs_1d_2d.py --X path/to/X.tsv --y path/to/y.tsv --cohort-name my_cohort [--out-dir results]
"""

import argparse
import os
import subprocess
import sys
import tempfile

import pandas as pd

# Same as main pipeline: 3 runs per cohort, then aggregate with mean
N_RUNS_PER_COHORT = 3

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)


def _find_run_mdfs_r():
    r_path = os.path.join(REPO_ROOT, "run_mdfs.R")
    if os.path.isfile(r_path):
        return r_path
    if os.path.isfile("run_mdfs.R"):
        return os.path.abspath("run_mdfs.R")
    return None


def run_mdfs_r(
    X_path: str,
    y_path: str,
    out_2d: str,
    out_unique: str,
    out_1d: str,
    seed: int = None,
    no_bonferroni: bool = False,
) -> None:
    """Call run_mdfs.R. If seed is set, pass it as 6th argument so each run gets different MDFS results.
    If no_bonferroni is True, pass 'no_bonferroni' as 7th argument to skip the strict
    pair-level Bonferroni threshold (used for exploratory synergy enumeration in apply_framework,
    while the CRC pipeline keeps the threshold for classifier feature selection)."""
    rscript = _find_run_mdfs_r()
    if not rscript:
        raise FileNotFoundError("run_mdfs.R not found in repo root or current directory")
    cmd = ["Rscript", rscript, X_path, y_path, out_2d, out_unique, out_1d]
    # seed must be arg 6, no_bonferroni must be arg 7
    cmd.append(str(seed) if seed is not None else "12")
    if no_bonferroni:
        cmd.append("no_bonferroni")
    result = subprocess.run(
        cmd,
        check=False,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"run_mdfs.R exited with code {result.returncode}.\n"
            f"STDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )


def _one_mdfs_run(X: pd.DataFrame, y: pd.Series, tmp: str, cohort_name: str, run_idx: int,
                   no_bonferroni: bool = False):
    """Run MDFS once with seed (12 + run_idx); return (ig_info_df, tuples_df, ig_1d_df)."""
    train_X = os.path.join(tmp, "X.tsv")
    train_y = os.path.join(tmp, "y.tsv")
    out_2d = os.path.join(tmp, "mdfs_2d.tsv")
    out_unique = os.path.join(tmp, "mdfs_unique.txt")
    out_1d = os.path.join(tmp, "mdfs_1d_ig.tsv")
    X.to_csv(train_X, sep="\t")
    y.rename("label").to_csv(train_y, sep="\t", header=True)
    run_seed = 12 + run_idx
    run_mdfs_r(train_X, train_y, out_2d, out_unique, out_1d, seed=run_seed,
               no_bonferroni=no_bonferroni)  # False by default: Bonferroni enabled

    # Read 2D tuples (new format: base, contributing, base_IG_1D, contributing_IG_1D, IG_2D_added, total_IG)
    tuples_df = pd.read_csv(out_2d, sep="\t")
    if tuples_df.empty:
        # No significant 2D features for this run
        ig_1d_df = pd.read_csv(out_1d, sep="\t")
        return pd.DataFrame(), tuples_df, ig_1d_df

    tuples_df["base"] = tuples_df["base"].astype(str).str.replace(" ", "_")
    tuples_df["contributing"] = tuples_df["contributing"].astype(str).str.replace(" ", "_")

    # Read 1D IGs (new format: Feature, IG, p_value, adjusted_p_value)
    ig_1d_df = pd.read_csv(out_1d, sep="\t")
    ig_1d_df["Feature"] = ig_1d_df["Feature"].astype(str).str.replace(" ", "_")

    # Build ig_info rows from tuples
    rows = []
    for _, r in tuples_df.iterrows():
        rows.append({
            "cohort": cohort_name,
            "run": run_idx,
            "base": r["base"],
            "contributing": r["contributing"],
            "base_IG_1D": float(r["base_IG_1D"]),
            "contributing_IG_1D": float(r["contributing_IG_1D"]),
            "IG_2D_added": float(r["IG_2D_added"]),
            "total_IG": float(r["total_IG"]),
        })
    return pd.DataFrame(rows), tuples_df, ig_1d_df


def run_mdfs_1d_2d(
    X: pd.DataFrame,
    y: pd.Series,
    cohort_name: str,
    out_dir: str,
    no_bonferroni: bool = False,
) -> None:
    """
    Run MDFS 1D and 2D on (X, y) N_RUNS_PER_COHORT times, aggregate IGs with mean,
    write outputs under out_dir.
    X: samples x features, y: binary labels (0=disease, 1=control/healthy).

    no_bonferroni: If False (default), apply the Bonferroni-corrected pair-level
        IG threshold, matching the CRC discovery pipeline (framework.ipynb).
    """
    os.makedirs(out_dir, exist_ok=True)
    safe_name = cohort_name.replace(os.sep, "_")

    with tempfile.TemporaryDirectory(prefix="mdfs_apply_") as tmp:
        all_ig = []
        tuples_df = None
        ig_1d_df = None
        for run_idx in range(N_RUNS_PER_COHORT):
            print(f"  Run {run_idx + 1}/{N_RUNS_PER_COHORT}...")
            ig_info, tuples_df, ig_1d_df = _one_mdfs_run(
                X, y, tmp, cohort_name, run_idx, no_bonferroni=no_bonferroni)
            if not ig_info.empty:
                all_ig.append(ig_info)

        # Write last-run 2D tuples and 1D IGs
        tuples_df.to_csv(os.path.join(out_dir, f"mdfs_2d_{safe_name}.tsv"), index=False, sep="\t")
        ig_1d_df.to_csv(os.path.join(out_dir, f"mdfs_1d_ig_{safe_name}.tsv"), index=False, sep="\t")

        if all_ig:
            ig_info = pd.concat(all_ig, ignore_index=True)
            ig_info.to_csv(os.path.join(out_dir, f"mdfs_ig_info_{safe_name}.tsv"), index=False, sep="\t")

            # Aggregated summary by (base, contributing): mean across runs
            summary = ig_info.groupby(["base", "contributing"]).agg(
                base_IG_1D_mean=("base_IG_1D", "mean"),
                contributing_IG_1D_mean=("contributing_IG_1D", "mean"),
                IG_2D_added_mean=("IG_2D_added", "mean"),
                IG_2D_added_std=("IG_2D_added", "std"),
                total_IG_mean=("total_IG", "mean"),
                n_runs=("IG_2D_added", "count"),
            ).reset_index()
            # IG gain: how much does total_IG exceed base_IG_1D?
            # This is simply IG_2D_added_mean (by definition), but we include it explicitly
            summary["ig_gain"] = summary["IG_2D_added_mean"]
            denom = summary["base_IG_1D_mean"].replace(0, float("nan"))
            summary["ig_gain_pct"] = (summary["ig_gain"] / denom) * 100

            summary.to_csv(os.path.join(out_dir, f"ig_gain_summary_{safe_name}.tsv"), index=False, sep="\t")

            print(
                f"  Wrote mdfs_2d_{safe_name}.tsv, mdfs_1d_ig_{safe_name}.tsv, "
                f"mdfs_ig_info_{safe_name}.tsv ({len(ig_info)} rows, {N_RUNS_PER_COHORT} runs), "
                f"ig_gain_summary_{safe_name}.tsv ({len(summary)} base-contributing pairs)"
            )
        else:
            print(f"  No significant 2D features across {N_RUNS_PER_COHORT} runs.")
            # Write empty ig_info and summary
            pd.DataFrame().to_csv(os.path.join(out_dir, f"mdfs_ig_info_{safe_name}.tsv"), index=False, sep="\t")
            pd.DataFrame().to_csv(os.path.join(out_dir, f"ig_gain_summary_{safe_name}.tsv"), index=False, sep="\t")
    return


def main():
    parser = argparse.ArgumentParser(
        description="Run MDFS 1D and 2D on cohort data and report IGs for pairs."
    )
    parser.add_argument(
        "--data-dir",
        default="data",
        help="Directory containing cohort data (merged or per-cohort subdirs)",
    )
    parser.add_argument(
        "--cohort",
        help="Cohort identifier (e.g. VogtmannE_2016). Ignored if --all-cohorts.",
    )
    parser.add_argument(
        "--all-cohorts",
        action="store_true",
        help="Run MDFS for every cohort found in data-dir",
    )
    parser.add_argument(
        "--X",
        help="Path to feature matrix TSV (samples x features). Use with --y and --cohort-name.",
    )
    parser.add_argument(
        "--y",
        help="Path to labels TSV (sample_id, label). Use with --X and --cohort-name.",
    )
    parser.add_argument(
        "--cohort-name",
        help="Cohort name for output files when using --X and --y",
    )
    parser.add_argument(
        "--out-dir",
        default="results",
        help="Directory for MDFS outputs (default: results)",
    )
    args = parser.parse_args()

    def _resolve(path):
        if os.path.isabs(path):
            return path
        return os.path.abspath(os.path.join(os.getcwd(), path))

    data_dir = _resolve(args.data_dir)
    out_dir = _resolve(args.out_dir)

    if args.X and args.y and args.cohort_name:
        X = pd.read_csv(args.X, sep="\t", index_col=0)
        y = pd.read_csv(args.y, sep="\t", index_col=0).squeeze()
        y = y.astype(int)
        common = X.index.intersection(y.index)
        X = X.loc[common]
        y = y.loc[common]
        print(f"Loaded {len(X)} samples, {X.shape[1]} features for cohort {args.cohort_name}")
        run_mdfs_1d_2d(X, y, args.cohort_name, out_dir)
        return

    try:
        from load_cohort_data import get_cohorts_with_data, load_combined_matrix_for_cohort
    except ImportError:
        sys.path.insert(0, SCRIPT_DIR)
        from load_cohort_data import get_cohorts_with_data, load_combined_matrix_for_cohort

    if not os.path.isdir(data_dir):
        print(f"Error: data directory not found: {data_dir}", file=sys.stderr)
        sys.exit(1)

    cohorts = get_cohorts_with_data(data_dir)
    if not cohorts:
        print("No cohorts found in data dir. Check metadata and abundance/pathway files.", file=sys.stderr)
        sys.exit(1)

    to_run = []
    if args.all_cohorts:
        to_run = list(cohorts)
    elif args.cohort:
        if args.cohort in cohorts:
            to_run = [args.cohort]
        else:
            print(f"Cohort '{args.cohort}' not found. Available: {list(cohorts)[:10]}...", file=sys.stderr)
            sys.exit(1)
    else:
        print("Specify --cohort <id>, --all-cohorts, or --X/--y/--cohort-name.", file=sys.stderr)
        sys.exit(1)

    failed = []
    for c in to_run:
        print(f"\n--- Cohort: {c} ---")
        try:
            X, y = load_combined_matrix_for_cohort(data_dir, c)
            if X is None or y is None or len(X) < 10:
                print(f"  Skip (no data or too few samples)")
                continue
            run_mdfs_1d_2d(X, y, c, out_dir)
        except Exception as e:
            print(f"  ERROR processing {c}: {e}", file=sys.stderr)
            failed.append(c)
            continue
    print("\nDone.")
    if failed:
        print(f"Failed cohorts ({len(failed)}): {failed}", file=sys.stderr)


if __name__ == "__main__":
    main()
