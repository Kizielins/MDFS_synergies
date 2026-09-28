#!/usr/bin/env python3
"""
Build synthetic features from MDFS synergistic pairs and evaluate them with a
Random Forest on a single cohort.

For each synergistic pair (f1, f2) two synthetic features are created:
  - LR_f1__f2 : log-ratio          log((f1 + eps) / (f2 + eps))
  - GM_f1__f2 : geometric mean     sqrt((f1 + eps) * (f2 + eps))

Pair selection (per MDFS run, repeated N_RUNS times with different seeds):
  1. Run MDFS 1D + 2D (run_mdfs.R) and keep significant, truly synergistic pairs.
  2. Score each unordered pair by the larger total_IG of its two directions.
  3. Keep the top --top-fraction of pairs (at least one).
  Consensus pairs = pairs kept in every run.

Evaluation (stratified k-fold cross-validation; MDFS is run on the training
folds only, so the test fold never informs pair selection). Three feature sets
are compared:
  - All Features               original features only
  - Synthetic Only             LR_/GM_ features only
  - All Features + Synthetic   both
For each, the number of features k is tuned on an internal validation split,
the top-k features by RF importance are selected, and a final RF is trained.

Usage:
    python synthetic_features.py --X samples_features.tsv --y labels.tsv --out-dir out
    python synthetic_features.py --X samples_features.tsv --y labels.tsv --out-dir out --no-cv

Input formats: as for compute_synergies.py
    X: tab-separated, first column = sample IDs, remaining columns = feature values
    y: tab-separated, first column = sample IDs, second column = labels (0=disease, 1=control)

Outputs (in --out-dir):
    synergy_pairs.tsv        consensus pairs selected on the full dataset
    synthetic_features.tsv   original features + LR_/GM_ features (full dataset)
    cv_performance.tsv       per-fold AUC / accuracy for each feature set (skipped with --no-cv)
    selected_features.tsv    how often each feature was selected across folds (skipped with --no-cv)
"""

import argparse
import os
import sys
import tempfile
import warnings
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split

from compute_synergies import N_RUNS_DEFAULT, _SEEDS, _one_mdfs_run

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# --- MODEL AND FEATURE SELECTION PARAMETERS ---
TOP_K_VALUES = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
MODEL_PARAMS = {
    'n_estimators': 1000, 'max_features': 'sqrt', 'min_samples_leaf': 1,
    'class_weight': 'balanced', 'random_state': 42, 'n_jobs': -1
}
SYNTHETIC_PREFIXES = ('LR_', 'GM_')
MIN_FEATURES = 10          # feature sets smaller than this are not evaluated
TOP_FRACTION_DEFAULT = 0.01
N_FOLDS_DEFAULT = 5
SCENARIOS = ['All Features', 'Synthetic Only', 'All Features + Synthetic']


def generate_synthetic_features(X, feature_pairs, epsilon=1e-9):
    """Return a DataFrame of LR_ and GM_ features for each (feature1, feature2) pair."""
    synthetic_X = pd.DataFrame(index=X.index)
    for pair in feature_pairs:
        f1, f2 = sorted(pair)
        if f1 in X.columns and f2 in X.columns:
            # Log-ratio (division)
            synthetic_X[f"LR_{f1}__{f2}"] = np.log((X[f1] + epsilon) / (X[f2] + epsilon))
            # Geometric mean
            synthetic_X[f"GM_{f1}__{f2}"] = np.sqrt((X[f1] + epsilon) * (X[f2] + epsilon))
    return synthetic_X


def find_best_k(X_train, y_train, k_values, model_params):
    """Choose the number of top-importance features k using an internal 75/25 split."""
    if X_train.shape[1] < min(k_values):
        return X_train.shape[1]

    X_train_internal, X_val_internal, y_train_internal, y_val_internal = train_test_split(
        X_train, y_train, test_size=0.25, random_state=42, stratify=y_train)

    ranking_model = RandomForestClassifier(**model_params).fit(X_train_internal, y_train_internal)
    sorted_features = pd.Series(ranking_model.feature_importances_,
                                index=X_train_internal.columns).sort_values(ascending=False).index

    best_auc, best_k = -1, k_values[0]
    for k in k_values:
        k = min(k, len(sorted_features))
        top_k_features = sorted_features[:k]
        model = RandomForestClassifier(**model_params).fit(X_train_internal[top_k_features], y_train_internal)
        y_pred_proba = model.predict_proba(X_val_internal.reindex(columns=top_k_features, fill_value=0))[:, 1]
        try:
            auc = roc_auc_score(y_val_internal, y_pred_proba)
            if auc > best_auc:
                best_auc, best_k = auc, k
        except ValueError:
            continue
    return best_k


def train_and_evaluate(X_train, y_train, X_test, y_test, k_values=None, model_params=None):
    """Tune k, select top-k features by RF importance, train a final RF and score it on the test set."""
    if k_values is None:
        k_values = TOP_K_VALUES
    if model_params is None:
        model_params = MODEL_PARAMS

    if X_train.empty or X_train.shape[1] < MIN_FEATURES:
        return {'auc': np.nan, 'accuracy': np.nan, 'best_k': 0, 'selected_features': []}

    best_k = find_best_k(X_train, y_train, k_values, model_params)

    ranking_model = RandomForestClassifier(**model_params).fit(X_train, y_train)
    sorted_features = pd.Series(ranking_model.feature_importances_, index=X_train.columns).sort_values(ascending=False)
    top_features = sorted_features.nlargest(best_k).index.tolist()

    final_model = RandomForestClassifier(**model_params).fit(X_train[top_features], y_train)
    X_test_aligned = X_test.reindex(columns=top_features, fill_value=0)
    y_pred_proba = final_model.predict_proba(X_test_aligned)[:, 1]
    y_pred = final_model.predict(X_test_aligned)

    return {
        'auc': roc_auc_score(y_test, y_pred_proba),
        'accuracy': accuracy_score(y_test, y_pred),
        'best_k': best_k,
        'selected_features': top_features,
    }


def select_synergy_pairs(X, y, n_runs=N_RUNS_DEFAULT, top_fraction=TOP_FRACTION_DEFAULT):
    """
    Run MDFS n_runs times, keep the top `top_fraction` of pairs per run (by total_IG)
    and return the consensus pairs (kept in every run) as a DataFrame with columns
    feature_1, feature_2, total_IG_mean, sorted by total_IG_mean descending.
    """
    kept_per_run = []
    with tempfile.TemporaryDirectory(prefix="mdfs_syn_") as tmp:
        for i, seed in enumerate(_SEEDS[:n_runs]):
            print(f"    MDFS run {i + 1}/{n_runs} (seed={seed})...", file=sys.stderr)
            pairs = _one_mdfs_run(X, y, tmp, seed)
            score = {}
            for r in pairs.itertuples(index=False):
                key = tuple(sorted((r.base, r.contributing)))
                score[key] = max(score.get(key, 0.0), r.total_IG)
            ranked = sorted(score.items(), key=lambda kv: -kv[1])
            n_keep = max(1, int(len(ranked) * top_fraction)) if ranked else 0
            kept_per_run.append(dict(ranked[:n_keep]))

    counts = Counter(p for run in kept_per_run for p in run)
    consensus = [p for p, c in counts.items() if c == n_runs]
    rows = [(a, b, np.mean([run[(a, b)] for run in kept_per_run])) for a, b in consensus]
    out = pd.DataFrame(rows, columns=['feature_1', 'feature_2', 'total_IG_mean'])
    return out.sort_values('total_IG_mean', ascending=False).reset_index(drop=True)


def cross_validate(X, y, n_folds, n_runs, top_fraction, seed=42):
    """Stratified k-fold CV of the three feature sets; pairs are selected on training folds only."""
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    perf_rows = []
    selected = Counter()
    for fold, (tr, te) in enumerate(skf.split(X, y), start=1):
        print(f"  Fold {fold}/{n_folds}", file=sys.stderr)
        X_tr, X_te, y_tr, y_te = X.iloc[tr], X.iloc[te], y.iloc[tr], y.iloc[te]

        pairs_df = select_synergy_pairs(X_tr, y_tr, n_runs, top_fraction)
        pairs = list(zip(pairs_df['feature_1'], pairs_df['feature_2']))
        syn_tr = generate_synthetic_features(X_tr, pairs)
        syn_te = generate_synthetic_features(X_te, pairs)

        scenarios = {
            'All Features': (X_tr, X_te),
            'Synthetic Only': (syn_tr, syn_te),
            'All Features + Synthetic': (pd.concat([X_tr, syn_tr], axis=1), pd.concat([X_te, syn_te], axis=1)),
        }
        for name, (a, b) in scenarios.items():
            res = train_and_evaluate(a, y_tr, b, y_te)
            perf_rows.append({'fold': fold, 'n_test': len(te), 'n_pairs': len(pairs), 'feature_set': name,
                              'n_features': a.shape[1], 'best_k': res['best_k'],
                              'auc': res['auc'], 'accuracy': res['accuracy']})
            if name == 'All Features + Synthetic':
                selected.update(res['selected_features'])

    perf = pd.DataFrame(perf_rows)
    sel = pd.DataFrame(selected.items(), columns=['feature', 'n_folds_selected'])
    sel['synthetic'] = sel['feature'].str.startswith(SYNTHETIC_PREFIXES)
    sel = sel.sort_values(['n_folds_selected', 'feature'], ascending=[False, True]).reset_index(drop=True)
    return perf, sel


def main():
    parser = argparse.ArgumentParser(
        description="Build LR_/GM_ synthetic features from MDFS synergistic pairs and "
                    "evaluate them with a Random Forest (single cohort, k-fold CV).",
    )
    parser.add_argument("--X", required=True,
                        help="Feature matrix TSV (samples x features; first column = sample IDs)")
    parser.add_argument("--y", required=True,
                        help="Labels TSV (two columns: sample_id, label; 0 = disease, 1 = control/healthy)")
    parser.add_argument("--out-dir", required=True, help="Directory for output tables")
    parser.add_argument("--n-runs", type=int, default=N_RUNS_DEFAULT,
                        help=f"MDFS runs per pair selection; consensus = pairs kept in all runs "
                             f"(default: {N_RUNS_DEFAULT}, max: {len(_SEEDS)})")
    parser.add_argument("--top-fraction", type=float, default=TOP_FRACTION_DEFAULT,
                        help=f"Fraction of pairs (by total_IG) kept per MDFS run, at least one "
                             f"(default: {TOP_FRACTION_DEFAULT})")
    parser.add_argument("--n-folds", type=int, default=N_FOLDS_DEFAULT,
                        help=f"Number of stratified cross-validation folds (default: {N_FOLDS_DEFAULT})")
    parser.add_argument("--no-cv", action="store_true",
                        help="Skip cross-validation; only select pairs and write synthetic features")
    args = parser.parse_args()

    if not 1 <= args.n_runs <= len(_SEEDS):
        parser.error(f"--n-runs must be between 1 and {len(_SEEDS)}")

    X = pd.read_csv(args.X, sep="\t", index_col=0)
    y = pd.read_csv(args.y, sep="\t", index_col=0).squeeze().astype(int)
    # MDFS pair names use underscores; keep X column names consistent with them
    X.columns = X.columns.astype(str).str.replace(" ", "_")

    common = X.index.intersection(y.index)
    if len(common) == 0:
        print("Error: no shared sample IDs between X and y.", file=sys.stderr)
        sys.exit(1)
    X, y = X.loc[common], y.loc[common]
    print(f"Loaded {len(X)} samples, {X.shape[1]} features.", file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)

    if not args.no_cv:
        print(f"Cross-validation ({args.n_folds} folds)...", file=sys.stderr)
        perf, sel = cross_validate(X, y, args.n_folds, args.n_runs, args.top_fraction)
        perf.to_csv(os.path.join(args.out_dir, "cv_performance.tsv"), sep="\t", index=False, float_format="%.6f")
        sel.to_csv(os.path.join(args.out_dir, "selected_features.tsv"), sep="\t", index=False)

        summary = perf.groupby('feature_set', sort=False)[['auc', 'accuracy']].agg(['mean', 'std'])
        print("\nCross-validated performance (mean ± sd across folds):", file=sys.stderr)
        for name in SCENARIOS:
            auc_m, auc_s, acc_m, acc_s = summary.loc[name]
            print(f"  {name:26s}  AUC = {auc_m:.3f} ± {auc_s:.3f}   accuracy = {acc_m:.3f} ± {acc_s:.3f}",
                  file=sys.stderr)

    print("\nSelecting synergy pairs on the full dataset...", file=sys.stderr)
    pairs_df = select_synergy_pairs(X, y, args.n_runs, args.top_fraction)
    pairs_df.to_csv(os.path.join(args.out_dir, "synergy_pairs.tsv"), sep="\t", index=False, float_format="%.6f")
    pairs = list(zip(pairs_df['feature_1'], pairs_df['feature_2']))
    X_out = pd.concat([X, generate_synthetic_features(X, pairs)], axis=1)
    X_out.to_csv(os.path.join(args.out_dir, "synthetic_features.tsv"), sep="\t")
    print(f"Wrote {len(pairs)} consensus pairs and {X_out.shape[1]} features to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
