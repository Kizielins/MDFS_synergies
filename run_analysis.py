#!/usr/bin/env python3
"""
MDFS synergy analysis of a single cohort.

Without arguments the bundled test data (test_files/) is analysed:
    python run_analysis.py
Your own data:
    python run_analysis.py --X my_features.tsv --y my_labels.tsv --out-dir my_results

Steps
  1. MDFS on the full dataset (--n-mdfs-runs runs): synergy table of all pairs
     passing the MDFS filters, and the consensus pairs (top --top-fraction of
     pairs in >= --min-runs runs) with their LR_/GM_ synthetic features.
  2. Model comparison on --n-splits stratified train/test splits. In each split
     MDFS is run --n-mdfs-runs times on the training part only, consensus pairs
     are turned into LR_/GM_ features, and two Random Forests are trained:
       Baseline RF  original features
       MDFS RF      original features + LR_/GM_ features of the consensus pairs
     Both are scored (AUC, accuracy) on the held-out test part.
  3. ROC curves of the two models and summary tables.

Input formats
  X: tab-separated, first column = sample IDs, other columns = feature values
  y: tab-separated with a header, columns sample_id and label (0 = disease, 1 = control)
"""

import argparse
import os
import sys
import warnings
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
from sklearn.model_selection import StratifiedShuffleSplit

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from mdfs_pairs import consensus_pairs, run_mdfs_repeated, synergy_table  # noqa: E402
from rf_model import SYNTHETIC_PREFIXES, generate_synthetic_features, train_and_evaluate  # noqa: E402

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_X = os.path.join(SCRIPT_DIR, "test_files", "X_test.tsv")
DEFAULT_Y = os.path.join(SCRIPT_DIR, "test_files", "y_test.tsv")

MODELS = ["Baseline RF", "MDFS RF"]
MODEL_COLORS = {"Baseline RF": "#eb6834", "MDFS RF": "#2a78d6"}


def read_tsv(path, what):
    if not os.path.isfile(path):
        sys.exit(f"Error: {what} file not found: {path}")
    df = pd.read_csv(path, sep="\t", index_col=0)
    if df.shape[1] == 0:
        sys.exit(f"Error: {what} file {path} has a single column. Is it tab-separated?")
    df.index = df.index.astype(str)
    if df.index.duplicated().any():
        sys.exit(f"Error: {what} file has duplicate sample IDs, e.g. {df.index[df.index.duplicated()][0]}.")
    return df


def load_data(x_path, y_path):
    X = read_tsv(x_path, "feature matrix (--X)")
    y_df = read_tsv(y_path, "labels (--y)")
    if y_df.shape[1] != 1:
        sys.exit(f"Error: the labels file must have two columns (sample_id, label); found {y_df.shape[1] + 1}.")
    if str(y_df.columns[0]).strip() in ("0", "1", "0.0", "1.0"):
        sys.exit("Error: the labels file seems to have no header row; the first line must be 'sample_id<TAB>label'.")

    # MDFS pair names use underscores; keep X column names consistent with them
    X.columns = X.columns.astype(str).str.replace(" ", "_")
    if X.columns.duplicated().any():
        sys.exit(f"Error: duplicate feature names (after replacing spaces with underscores), "
                 f"e.g. {X.columns[X.columns.duplicated()][0]}.")
    non_numeric = [c for c in X.columns if not pd.api.types.is_numeric_dtype(X[c])]
    if non_numeric:
        sys.exit(f"Error: non-numeric values in feature(s) {non_numeric[:3]}.")
    if X.isna().any().any():
        sys.exit("Error: the feature matrix contains missing values.")
    if (X < 0).any().any():
        sys.exit("Error: the feature matrix contains negative values. Log-ratio and geometric-mean features "
                 "need non-negative data such as (relative) abundances, not CLR- or log-transformed values.")

    labels = pd.to_numeric(y_df.iloc[:, 0], errors="coerce")
    bad = labels.isna() | ~labels.isin([0, 1])
    if bad.any():
        examples = sorted(set(y_df.iloc[:, 0][bad].astype(str)))[:5]
        sys.exit(f"Error: labels must be 0 (disease) or 1 (control); found {examples}.")
    y = labels.astype(int)

    common = X.index.intersection(y.index)
    if len(common) == 0:
        sys.exit("Error: no shared sample IDs between X and y.")
    only_x, only_y = len(X.index.difference(y.index)), len(y.index.difference(X.index))
    if only_x or only_y:
        print(f"Note: {only_x} samples only in X and {only_y} only in y are ignored.", file=sys.stderr)
    X, y = X.loc[common], y.loc[common]
    if y.nunique() < 2:
        sys.exit("Error: both classes (0 and 1) must be present.")
    return X, y


def display_path(path):
    """Path for messages: relative to the repository for bundled files, as given otherwise."""
    full = os.path.abspath(path)
    return os.path.relpath(full, SCRIPT_DIR) if full.startswith(SCRIPT_DIR + os.sep) else path


def evaluate_split(X, y, train_idx, test_idx, args, split):
    """Consensus MDFS pairs on the training part, then baseline and MDFS RF on the test part."""
    X_tr, X_te, y_tr, y_te = X.iloc[train_idx], X.iloc[test_idx], y.iloc[train_idx], y.iloc[test_idx]

    runs = run_mdfs_repeated(X_tr, y_tr, args.n_mdfs_runs, args.seed, label=f"split {split}: ")
    pairs_df = consensus_pairs(runs, args.top_fraction, args.min_runs)
    pairs = list(zip(pairs_df["feature_1"], pairs_df["feature_2"]))

    feature_sets = {
        "Baseline RF": (X_tr, X_te),
        "MDFS RF": (pd.concat([X_tr, generate_synthetic_features(X_tr, pairs)], axis=1),
                    pd.concat([X_te, generate_synthetic_features(X_te, pairs)], axis=1)),
    }
    results = {}
    for model, (a, b) in feature_sets.items():
        res = train_and_evaluate(a, y_tr, b, y_te)
        res["n_features"] = a.shape[1]
        results[model] = res
    return pairs, results, y_te.to_numpy()


def plot_roc(roc_data, auc_table, out_base):
    """Mean ROC curve (± sd band) per model across splits, and paired per-split AUCs."""
    grid = np.linspace(0, 1, 201)
    fig, (ax_roc, ax_auc) = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw={"width_ratios": [1.3, 1]})

    for model in MODELS:
        color = MODEL_COLORS[model]
        tprs = []
        for fpr, tpr in roc_data[model]:
            ax_roc.plot(fpr, tpr, color=color, alpha=0.15, lw=1)
            interp = np.interp(grid, fpr, tpr)
            interp[0] = 0.0
            tprs.append(interp)
        mean_tpr, sd_tpr = np.mean(tprs, axis=0), np.std(tprs, axis=0)
        mean_tpr[-1] = 1.0
        aucs = auc_table.loc[auc_table["model"] == model, "auc"]
        ax_roc.plot(grid, mean_tpr, color=color, lw=2,
                    label=f"{model}  AUC = {aucs.mean():.3f} ± {aucs.std():.3f}")
        ax_roc.fill_between(grid, np.clip(mean_tpr - sd_tpr, 0, 1), np.clip(mean_tpr + sd_tpr, 0, 1),
                            color=color, alpha=0.15, lw=0)
    ax_roc.plot([0, 1], [0, 1], ls="--", color="#8a8984", lw=1)
    ax_roc.set(xlim=(0, 1), ylim=(0, 1.01), xlabel="False positive rate", ylabel="True positive rate",
               title=f"ROC curves (mean ± sd over {len(roc_data[MODELS[0]])} test splits)")
    ax_roc.legend(loc="lower right", frameon=False)

    wide = auc_table.pivot(index="split", columns="model", values="auc")[MODELS]
    for _, row in wide.iterrows():
        ax_auc.plot([0, 1], row.to_numpy(), color="#b5b4ae", lw=1, zorder=1)
    for i, model in enumerate(MODELS):
        ax_auc.scatter(np.full(len(wide), i), wide[model], s=40, color=MODEL_COLORS[model],
                       edgecolor="white", linewidth=1.5, zorder=2)
    ax_auc.set(xticks=[0, 1], xticklabels=MODELS, xlim=(-0.4, 1.4), ylabel="Test AUC",
               title="Test AUC per split")

    for ax in (ax_roc, ax_auc):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", dpi=200)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="MDFS synergy analysis of a single cohort: synergy table, consensus pairs, "
                    "and Random Forest comparison (baseline vs MDFS synthetic features). "
                    "Runs on the bundled test data when --X/--y are not given.",
    )
    parser.add_argument("--X", default=DEFAULT_X,
                        help="Feature matrix TSV: samples x features, first column = sample IDs "
                             "(default: test_files/X_test.tsv)")
    parser.add_argument("--y", default=DEFAULT_Y,
                        help="Labels TSV: sample_id, label (0 = disease, 1 = control) "
                             "(default: test_files/y_test.tsv)")
    parser.add_argument("--out-dir", default="results", help="Output directory (default: results)")
    parser.add_argument("--n-mdfs-runs", type=int, default=10, help="MDFS runs per selection (default: 10)")
    parser.add_argument("--min-runs", type=int, default=6,
                        help="A pair is a consensus pair if selected in at least this many runs (default: 6)")
    parser.add_argument("--top-fraction", type=float, default=0.02,
                        help="Fraction of pairs (by total_IG) selected per MDFS run, at least one (default: 0.02)")
    parser.add_argument("--n-splits", type=int, default=10, help="Stratified train/test splits (default: 10)")
    parser.add_argument("--test-size", type=float, default=0.3, help="Test fraction per split (default: 0.3)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Seed of the train/test splits; MDFS run i uses seed + i (default: 42)")
    args = parser.parse_args()

    if not 1 <= args.min_runs <= args.n_mdfs_runs:
        parser.error("--min-runs must be between 1 and --n-mdfs-runs")
    if not 0 < args.top_fraction <= 1:
        parser.error("--top-fraction must be in (0, 1]")
    if not 0 < args.test_size < 1:
        parser.error("--test-size must be in (0, 1)")
    if args.n_splits < 2:
        parser.error("--n-splits must be at least 2")
    if args.n_mdfs_runs < 1:
        parser.error("--n-mdfs-runs must be at least 1")

    X, y = load_data(args.X, args.y)
    counts = y.value_counts()
    print(f"Loaded {len(X)} samples ({counts.get(0, 0)} disease, {counts.get(1, 0)} control), "
          f"{X.shape[1]} features from {display_path(args.X)}", file=sys.stderr)
    os.makedirs(args.out_dir, exist_ok=True)
    out = lambda name: os.path.join(args.out_dir, name)  # noqa: E731

    # ---- 1. MDFS on the full dataset ----
    print("\n[1/3] MDFS on the full dataset", file=sys.stderr)
    runs = run_mdfs_repeated(X, y, args.n_mdfs_runs, args.seed, label="full data: ")
    table = synergy_table(runs, args.top_fraction, args.min_runs)
    table.to_csv(out("mdfs_synergies.tsv"), sep="\t", index=False, float_format="%.6f")
    full_pairs = consensus_pairs(runs, args.top_fraction, args.min_runs)
    full_pairs.to_csv(out("consensus_pairs.tsv"), sep="\t", index=False, float_format="%.6f")
    pairs = list(zip(full_pairs["feature_1"], full_pairs["feature_2"]))
    pd.concat([X, generate_synthetic_features(X, pairs)], axis=1).to_csv(out("synthetic_features.tsv"), sep="\t")
    print(f"  {len(table)} significant synergistic pair directions; {len(full_pairs)} consensus pairs",
          file=sys.stderr)

    # ---- 2. Baseline vs MDFS RF on stratified splits ----
    print(f"\n[2/3] Baseline vs MDFS RF on {args.n_splits} stratified train/test splits", file=sys.stderr)
    splitter = StratifiedShuffleSplit(n_splits=args.n_splits, test_size=args.test_size, random_state=args.seed)
    auc_rows, roc_data = [], {m: [] for m in MODELS}
    split_pair_counts, selected_counts = Counter(), Counter()
    for split, (tr, te) in enumerate(splitter.split(X, y), start=1):
        split_pairs, results, y_te = evaluate_split(X, y, tr, te, args, split)
        split_pair_counts.update(split_pairs)
        selected_counts.update(results["MDFS RF"]["selected_features"])
        for model in MODELS:
            res = results[model]
            fpr, tpr, _ = roc_curve(y_te, res["y_score"])
            roc_data[model].append((fpr, tpr))
            auc_rows.append({"split": split, "model": model, "n_train": len(tr), "n_test": len(te),
                             "n_consensus_pairs": len(split_pairs), "n_features": res["n_features"],
                             "best_k": res["best_k"], "auc": res["auc"], "accuracy": res["accuracy"]})
        a, b = results["Baseline RF"]["auc"], results["MDFS RF"]["auc"]
        print(f"  split {split}: {len(split_pairs)} consensus pairs; AUC baseline {a:.3f}, MDFS {b:.3f}",
              file=sys.stderr)

    # ---- 3. Outputs ----
    print("\n[3/3] Writing outputs", file=sys.stderr)
    auc_table = pd.DataFrame(auc_rows)
    auc_table.to_csv(out("auc_per_split.tsv"), sep="\t", index=False, float_format="%.6f")

    split_pairs_df = pd.DataFrame([(a, b, n) for (a, b), n in split_pair_counts.items()],
                                  columns=["feature_1", "feature_2", "n_splits"])
    split_pairs_df.sort_values(["n_splits", "feature_1", "feature_2"], ascending=[False, True, True]) \
        .to_csv(out("split_consensus_pairs.tsv"), sep="\t", index=False)

    sel = pd.DataFrame(selected_counts.items(), columns=["feature", "n_splits_selected"])
    sel["synthetic"] = sel["feature"].str.startswith(SYNTHETIC_PREFIXES)
    sel.sort_values(["n_splits_selected", "feature"], ascending=[False, True]) \
        .to_csv(out("selected_features.tsv"), sep="\t", index=False)

    plot_roc(roc_data, auc_table, out("roc_comparison"))

    wide = auc_table.pivot(index="split", columns="model", values="auc")[MODELS]
    diff = wide["MDFS RF"] - wide["Baseline RF"]
    lines = [
        f"Input: {display_path(args.X)} ({len(X)} samples, {X.shape[1]} features); "
        f"labels: {display_path(args.y)}",
        f"MDFS: {args.n_mdfs_runs} runs (seeds {args.seed}-{args.seed + args.n_mdfs_runs - 1}), "
        f"top {args.top_fraction * 100:g}% of pairs per run, consensus = selected in >= {args.min_runs} runs",
        f"Full dataset: {len(table)} significant pair directions, {len(full_pairs)} consensus pairs",
        f"Evaluation: {args.n_splits} stratified splits, test size {args.test_size}",
        "",
    ]
    for model in MODELS:
        m = auc_table[auc_table["model"] == model]
        lines.append(f"{model:12s} AUC {m['auc'].mean():.3f} ± {m['auc'].std():.3f}   "
                     f"accuracy {m['accuracy'].mean():.3f} ± {m['accuracy'].std():.3f}")
    lines.append(f"AUC difference (MDFS - baseline): {diff.mean():+.3f} ± {diff.std():.3f}; "
                 f"MDFS RF higher in {(diff > 0).sum()}/{len(diff)} splits")
    summary = "\n".join(lines) + "\n"
    with open(out("summary.txt"), "w") as fh:
        fh.write(summary)
    print("\n" + summary + f"Results written to {args.out_dir}/", file=sys.stderr)


if __name__ == "__main__":
    main()
