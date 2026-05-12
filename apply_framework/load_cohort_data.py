#!/usr/bin/env python3
"""
Load cohort data from apply_framework data dir: combined taxa + functions matrix,
cohort IDs, and disease labels from metadata.

By default, CRC cohorts that were used to build the original pipeline are excluded
so that apply_framework runs on different cohorts. Use --include-crc-cohorts to include them.

Supports two layouts:
  1) Merged: data_dir/merged_abundances.txt (taxa), merged_pathways.txt (pathways),
     merged_metadata.tsv. Cohorts = unique study_name or source_cohort in metadata.
  2) Per-cohort: data_dir/taxa/<cohort>.tsv, data_dir/pathways/<cohort>.tsv,
     data_dir/metadata/<cohort>.tsv.

Disease labels: from 'study_condition' or 'disease'.
  - control / healthy → 1
  - All other conditions (CRC, Adenoma, or any other disease) → 0
Samples with missing or unrecognized labels are dropped.

Usage:
  python load_cohort_data.py --data-dir data [--mdfs-results-dir results] [--out-table synergy_by_cohort.tsv]
  python load_cohort_data.py --data-dir data --include-crc-cohorts   # to include pipeline CRC cohorts
"""

import argparse
import os
import sys

import pandas as pd


# Default column names for cohort and disease (can be overridden by metadata)
COHORT_COLUMNS = ["source_cohort", "study_name"]
LABEL_COLUMNS = ["study_condition", "disease"]
# Encode: control/healthy = 1, all other diseases = 0
# (CRC, Adenoma, or any other disease label → 0)
CONTROL_VALUES = {"control", "healthy"}

# CRC cohorts used to build the original pipeline — excluded so apply_framework uses different cohorts.
# Names with a prefix (e.g. Public_study__ZellerG_2014) are matched via the stem after the last "__".
EXCLUDED_CRC_COHORTS = {
    "FengQ_2015", "GuptaA_2019", "LiuNN_2022", "ThomasAM_2018b", "VogtmannE_2016",
    "WirbelJ_2018", "YachidaS_2019", "YangJ_2020", "YangY_2021", "YuJ_2015", "ZellerG_2014",
    "HanniganGD_2017", "BaxterNT_2016",
}


def _is_excluded_cohort(cohort_id: str) -> bool:
    """True if this cohort is one of the CRC cohorts used to build the pipeline."""
    c = str(cohort_id).strip()
    if c in EXCLUDED_CRC_COHORTS:
        return True
    # Match stem after last '__' (e.g. Public_study__ZellerG_2014 -> ZellerG_2014)
    if "__" in c:
        stem = c.split("__")[-1]
        if stem in EXCLUDED_CRC_COHORTS:
            return True
    return False


def _normalize_feature_names(series):
    return series.astype(str).str.replace(" ", "_", regex=False)


def get_cohorts_with_data(data_dir: str, exclude_crc_pipeline_cohorts: bool = True) -> list:
    """
    Discover cohort IDs from data_dir.
    Prefers merged_metadata.tsv (use study_name or source_cohort); else lists
    metadata files in metadata/ subdir (stem = cohort).
    By default excludes CRC cohorts that were used to build the original pipeline.
    """
    data_dir = os.path.abspath(data_dir)
    merged_meta = os.path.join(data_dir, "merged_metadata.tsv")
    if os.path.isfile(merged_meta):
        meta = pd.read_csv(merged_meta, sep="\t", low_memory=False)
        for col in COHORT_COLUMNS:
            if col in meta.columns:
                cohorts = meta[col].dropna().astype(str).unique().tolist()
                raw = sorted(set(c for c in cohorts if c and c != "nan"))
                if exclude_crc_pipeline_cohorts:
                    raw = [c for c in raw if not _is_excluded_cohort(c)]
                return raw
    meta_dir = os.path.join(data_dir, "metadata")
    if os.path.isdir(meta_dir):
        cohorts = []
        for f in os.listdir(meta_dir):
            if f.endswith(".tsv") or f.endswith(".csv"):
                cohorts.append(os.path.splitext(f)[0])
        if exclude_crc_pipeline_cohorts:
            cohorts = [c for c in cohorts if not _is_excluded_cohort(c)]
        return sorted(cohorts)
    return []


def _load_merged_abundance_pathways(data_dir: str):
    """Load merged taxa and pathways; return (taxa_df, pathways_df) or (None, None)."""
    taxa_path = os.path.join(data_dir, "merged_abundances.txt")
    path_path = os.path.join(data_dir, "merged_pathways.txt")
    if not os.path.isfile(taxa_path):
        return None, None
    taxa = pd.read_csv(taxa_path, sep="\t", index_col=0)
    # First column might be "Species" or feature names
    if taxa.index.name is None and len(taxa.columns) > 0:
        if taxa.columns[0].lower() in ("species", "feature", "pathway"):
            taxa = pd.read_csv(taxa_path, sep="\t", index_col=0)
    pathways = None
    if os.path.isfile(path_path):
        pathways = pd.read_csv(path_path, sep="\t", index_col=0)
        pathways.index = _normalize_feature_names(pd.Series(pathways.index))
    return taxa, pathways


def _encode_labels(labels: pd.Series) -> pd.Series:
    """Encode: control/healthy → 1, any other (CRC, Adenoma, etc.) → 0. Unrecognized/NaN → -1 (drop)."""
    out = []
    for v in labels.astype(str).str.strip().str.lower():
        if v in CONTROL_VALUES:
            out.append(1)
        elif v and v not in ("nan", ""):
            out.append(0)  # any disease
        else:
            out.append(-1)
    return pd.Series(out, index=labels.index)


def load_combined_matrix_for_cohort(data_dir: str, cohort_id: str):
    """
    Load one cohort's X (samples x features) and y (binary 0/1) with taxa + functions
    combined. Returns (X, y) or (None, None).
    """
    data_dir = os.path.abspath(data_dir)

    # --- Per-cohort files ---
    meta_path = os.path.join(data_dir, "metadata", f"{cohort_id}.tsv")
    if not os.path.isfile(meta_path):
        meta_path = os.path.join(data_dir, "metadata", f"{cohort_id}.csv")
    taxa_path = os.path.join(data_dir, "taxa", f"{cohort_id}.tsv")
    path_path = os.path.join(data_dir, "pathways", f"{cohort_id}.tsv")

    if os.path.isfile(meta_path) and os.path.isfile(taxa_path):
        return _load_per_cohort_files(data_dir, cohort_id, meta_path, taxa_path, path_path)

    # --- Merged layout ---
    merged_meta = os.path.join(data_dir, "merged_metadata.tsv")
    if not os.path.isfile(merged_meta):
        return None, None
    meta = pd.read_csv(merged_meta, sep="\t", low_memory=False)
    cohort_col = None
    for c in COHORT_COLUMNS:
        if c in meta.columns:
            cohort_col = c
            break
    if cohort_col is None:
        return None, None
    meta_cohort = meta[meta[cohort_col].astype(str) == str(cohort_id)]
    if meta_cohort.empty:
        return None, None

    label_col = None
    for c in LABEL_COLUMNS:
        if c in meta_cohort.columns:
            label_col = c
            break
    if label_col is None:
        return None, None

    sample_ids = set(meta_cohort["sample_id"].astype(str))
    taxa, pathways = _load_merged_abundance_pathways(data_dir)
    if taxa is None:
        return None, None

    # Align to samples present in metadata
    taxa_cols = [c for c in taxa.columns if str(c).strip() in sample_ids]
    if not taxa_cols:
        return None, None
    path_cols = []
    if pathways is not None and not pathways.empty:
        path_cols = [c for c in pathways.columns if str(c).strip() in sample_ids]
    common_samples = list(set(taxa_cols) & set(path_cols)) if path_cols else taxa_cols

    X_taxa = taxa[common_samples].T  # samples x features
    X_taxa.index = X_taxa.index.astype(str)
    if pathways is not None and not pathways.empty and common_samples:
        path_common = [c for c in common_samples if c in pathways.columns]
        if path_common:
            X_path = pathways[path_common].T
            X_path.index = X_path.index.astype(str)
            common_idx = X_taxa.index.intersection(X_path.index)
            X_taxa = X_taxa.loc[common_idx]
            X_path = X_path.reindex(common_idx, fill_value=0).fillna(0)
            X = pd.concat([X_taxa, X_path], axis=1)
        else:
            X = X_taxa.copy()
    else:
        X = X_taxa.copy()

    meta_cohort = meta_cohort.set_index("sample_id")
    meta_cohort.index = meta_cohort.index.astype(str)
    common = X.index.intersection(meta_cohort.index)
    X = X.loc[common]
    y_raw = meta_cohort.loc[common, label_col]
    y = _encode_labels(y_raw)
    # Drop non-binary (e.g. Adenoma)
    valid = y >= 0
    X = X[valid]
    y = y[valid]
    if y.nunique() < 2 or len(X) < 10:
        return None, None
    # Relative abundance
    X = X.div(X.sum(axis=1), axis=0).fillna(0)
    return X, y


def _load_per_cohort_files(data_dir, cohort_id, meta_path, taxa_path, path_path):
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    if meta.index.name is None and "sample_id" in meta.columns:
        meta = meta.set_index("sample_id")
    label_col = None
    for c in LABEL_COLUMNS:
        if c in meta.columns:
            label_col = c
            break
    if label_col is None:
        if "Study condition" in meta.columns:
            label_col = "Study condition"
        else:
            return None, None

    species = pd.read_csv(taxa_path, sep="\t", index_col=0)
    if "s__" in str(species.index[0]):
        drop = [i for i in species.index if "s__" not in str(i) or "t__" in str(i)]
        species = species.drop(drop, errors="ignore")
        species.index = [str(i).split("s__")[-1] for i in species.index]
    species = species.div(species.sum(axis=0), axis=1)
    X = species.T  # samples x features

    if os.path.isfile(path_path):
        paths = pd.read_csv(path_path, sep="\t", index_col=0)
        paths.index = _normalize_feature_names(pd.Series(paths.index))
        paths = paths.div(paths.sum(axis=0), axis=1)
        X_path = paths.T
        common_idx = X.index.intersection(X_path.index)
        X = X.loc[common_idx]
        X_path = X_path.reindex(common_idx, fill_value=0).fillna(0)
        X = pd.concat([X, X_path.reindex(columns=paths.index, fill_value=0)], axis=1)

    common = X.index.intersection(meta.index)
    X = X.loc[common]
    y_raw = meta.loc[common, label_col]
    y = _encode_labels(y_raw)
    valid = y >= 0
    X = X[valid]
    y = y[valid]
    if y.nunique() < 2 or len(X) < 10:
        return None, None
    X = X.div(X.sum(axis=1), axis=0).fillna(0)
    return X, y


def build_synergy_table(mdfs_results_dir: str, out_path: str = None) -> pd.DataFrame:
    """
    For each mdfs_ig_info_<cohort>.tsv in mdfs_results_dir, aggregate by
    (cohort, base, contributing) using mean across runs, then compute the IG gain:
    how much total_IG exceeds base_IG_1D (i.e., how much the contributing variable
    improves the base feature's information gain).

    New column format: cohort, run, base, contributing, base_IG_1D,
    contributing_IG_1D, IG_2D_added, total_IG.
    """
    if not os.path.isdir(mdfs_results_dir):
        return pd.DataFrame()
    rows = []
    for f in os.listdir(mdfs_results_dir):
        if not f.startswith("mdfs_ig_info_") or not f.endswith(".tsv"):
            continue
        cohort = f.replace("mdfs_ig_info_", "").replace(".tsv", "")
        path = os.path.join(mdfs_results_dir, f)
        if os.path.getsize(path) == 0:
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        if "cohort" not in df.columns:
            df["cohort"] = cohort
        rows.append(df)
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    # Aggregate by (cohort, base, contributing) with mean (handles multiple runs)
    group_cols = ["cohort", "base", "contributing"]
    agg_df = out.groupby(group_cols, as_index=False).agg(
        base_IG_1D_mean=("base_IG_1D", "mean"),
        contributing_IG_1D_mean=("contributing_IG_1D", "mean"),
        IG_2D_added_mean=("IG_2D_added", "mean"),
        total_IG_mean=("total_IG", "mean"),
        n_runs=("IG_2D_added", "count"),
    )
    # IG gain: total_IG vs base_IG_1D (how much does the contributing variable help?)
    agg_df["ig_gain"] = agg_df["IG_2D_added_mean"]
    denom = agg_df["base_IG_1D_mean"].replace(0, float("nan"))
    agg_df["ig_gain_pct"] = (agg_df["ig_gain"] / denom) * 100
    if out_path:
        agg_df.to_csv(out_path, index=False, sep="\t")
    return agg_df


def main():
    parser = argparse.ArgumentParser(
        description="Load cohort data (taxa+functions, labels) and build synergy table from MDFS results."
    )
    parser.add_argument("--data-dir", default="data", help="Data directory")
    parser.add_argument("--mdfs-results-dir", default="results", help="Directory with mdfs_ig_info_*.tsv")
    parser.add_argument("--out-table", default="synergy_by_cohort.tsv", help="Output table path")
    parser.add_argument("--list-cohorts", action="store_true", help="Only list discovered cohorts")
    parser.add_argument(
        "--include-crc-cohorts",
        action="store_true",
        help="Include CRC cohorts used to build the pipeline (default: exclude them)",
    )
    args = parser.parse_args()

    data_dir = os.path.join(os.path.dirname(__file__), args.data_dir)
    if not os.path.isdir(data_dir):
        data_dir = os.path.abspath(args.data_dir)
    if not os.path.isdir(data_dir):
        print("Data dir not found:", data_dir, file=sys.stderr)
        sys.exit(1)

    cohorts = get_cohorts_with_data(data_dir, exclude_crc_pipeline_cohorts=not args.include_crc_cohorts)
    print("Cohorts (CRC pipeline cohorts excluded by default):", len(cohorts))
    for c in cohorts[:20]:
        print(" ", c)
    if len(cohorts) > 20:
        print(" ... and", len(cohorts) - 20, "more")

    if args.list_cohorts:
        return

    res_dir = os.path.join(os.path.dirname(__file__), args.mdfs_results_dir)
    if not os.path.isdir(res_dir):
        res_dir = os.path.abspath(args.mdfs_results_dir)
    out_path = os.path.join(os.path.dirname(__file__), args.out_table)
    table = build_synergy_table(res_dir, out_path)
    print("IG gain table:", len(table), "base-contributing pairs")
    if not table.empty:
        print("Written to", out_path)


if __name__ == "__main__":
    main()
