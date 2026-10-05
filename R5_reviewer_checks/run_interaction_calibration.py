#!/usr/bin/env python3
"""Null-calibrated interaction information (Reviewer 2, comment 2; Reviewer 1, comment 7).

Interaction information of a pair: II = total_IG - IG_1D(feature 1) - IG_1D(feature 2)
(> 0: the pair carries more information about the label than its two members add up to).
The estimate is biased upwards, so its chance level is measured: for every LOCO fold and data type
every feature is shuffled independently within each class (univariate signal kept, all dependence
between features removed; N_PERM permutations, MDFS seed 42) and run_mdfs.R is applied. The 99th
percentile of II among the pairs that pass the pipeline's filters in these null runs is the fold's
threshold. A pair is called super-additive in a fold-run if its II exceeds that threshold.

Outputs (this folder):
  interaction_null_thresholds.tsv   per data type and fold: null 95th/99th percentile, share of null pairs with II > 0
  interaction_calibrated_summary.tsv per data type: consensus pairs, share super-additive (II > threshold in >= 6 of 10 runs)
  pair_interaction_<dt>.tsv         per unordered pair: mean II, fold-runs observed, % of them super-additive
  stable_pairs_interaction.tsv      the stable pairs (Results 2) with these values
  r5_numbers_interaction.tsv        numbers for manuscript/RESULTS_NUMBERS.md
and it adds the columns interaction_information_mean, pct_fold_runs_super_additive and super_additive to
manuscript/SuppTable2_top30_synergistic_pairs.tsv and manuscript/SuppTable4_stable_pairs.tsv.
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import training_set, run_mdfs, HERE, REPO  # noqa: E402

N_PERM, SEED, Q = 3, 42, 99
RES = os.path.join(REPO, "M2_loco_crc_pipeline", "results")
MS = os.path.join(REPO, "manuscript")


def null_fold(dt, fold, test_cohort):
    X, y = training_set(dt, test_cohort)
    vals = []
    for perm in range(1, N_PERM + 1):
        rng = np.random.default_rng(1000 * fold + perm)
        arr = X.values.copy()
        for cls in (0, 1):
            idx = np.where(y.values == cls)[0]
            for j in range(arr.shape[1]):
                arr[idx, j] = arr[rng.permutation(idx), j]
        pairs, *_ = run_mdfs(pd.DataFrame(arr, index=X.index, columns=X.columns), y, SEED)
        if len(pairs):
            pairs["key"] = [tuple(sorted((a, b))) for a, b in zip(pairs.base, pairs.contributing)]
            p = pairs.drop_duplicates("key")
            vals.append((p.total_IG - p.base_IG_1D - p.contributing_IG_1D).values)
    v = np.concatenate(vals) if vals else np.array([0.0])
    return {"data_type": dt, "fold": fold, "test_cohort": test_cohort, "null_pairs_per_run": int(len(v) / N_PERM),
            "null_share_II_positive": float((v > 0).mean()), "null_median_II": float(np.median(v)),
            "null_q95": float(np.percentile(v, 95)), "null_q99": float(np.percentile(v, Q)), "null_max": float(v.max())}


def main():
    thr_path = os.path.join(HERE, "interaction_null_thresholds.tsv")
    if os.path.exists(thr_path) and "--recompute" not in sys.argv:
        thr = pd.read_csv(thr_path, sep="\t")          # reuse the saved null thresholds (pass --recompute to rerun the permutations)
    else:
        tasks = []
        for dt in ("taxa", "functions", "combined"):
            f = pd.read_csv(os.path.join(RES, dt, f"mdfs_ig_info_{dt}.csv"), usecols=["fold", "test_cohort"]).drop_duplicates()
            tasks += [(dt, int(a), b) for a, b in zip(f.fold, f.test_cohort)]
        rows = []
        with ProcessPoolExecutor(8) as ex:
            futs = [ex.submit(null_fold, *t) for t in tasks]
            for i, fu in enumerate(as_completed(futs), 1):
                rows.append(fu.result()); print(f"[{i}/{len(tasks)}] {rows[-1]['data_type']} fold {rows[-1]['fold']}: null q99 = {rows[-1]['null_q99']:.2f}", flush=True)
        thr = pd.DataFrame(rows).sort_values(["data_type", "fold"]); thr.to_csv(os.path.join(HERE, "interaction_null_thresholds.tsv"), sep="\t", index=False)

    numbers, summ, per_pair = [], [], {}
    for dt in ("taxa", "functions", "combined"):
        d = pd.read_csv(os.path.join(RES, dt, f"mdfs_ig_info_{dt}.csv"))
        d["II"] = d.total_IG - d.base_IG_1D - d.contributing_IG_1D
        d["key"] = [tuple(sorted((a, b))) for a, b in zip(d.base, d.contributing)]
        u = d.drop_duplicates(["fold", "run", "key"]).merge(thr[thr.data_type == dt][["fold", "null_q99"]], on="fold")
        u["super"] = u.II > u.null_q99
        g = u.groupby(["fold", "key"]).agg(n_runs=("run", "nunique"), n_super=("super", "sum"), med_II=("II", "median")).reset_index()
        cons = g[g.n_runs >= 6]
        sa = (cons.n_super >= 6)
        summ.append({"data_type": dt, "consensus_pair_fold_entries": len(cons), "super_additive": int(sa.sum()), "share_super_additive": sa.mean(),
                     "share_median_II_negative": (cons.med_II < 0).mean(), "median_null_q99_bits": thr[thr.data_type == dt].null_q99.median(),
                     "null_share_II_positive": thr[thr.data_type == dt].null_share_II_positive.mean()})
        numbers.append({"section": "R5", "claim": f"{dt}: consensus pairs that are super-additive (interaction information above the within-class null 99th percentile in >= 6 of 10 runs)",
                        "value": f"{100 * sa.mean():.1f} % of {len(cons):,} pair-fold entries (null threshold median {thr[thr.data_type == dt].null_q99.median():.1f} bits; under the null {100 * thr[thr.data_type == dt].null_share_II_positive.mean():.0f} % of pairs have II > 0)",
                        "source": "R5 interaction_calibrated_summary.tsv"})
        pp = u.groupby("key").agg(interaction_information_mean=("II", "mean"), fold_runs_observed=("II", "size"), n_super=("super", "sum")).reset_index()
        pp["pct_fold_runs_super_additive"] = 100 * pp.n_super / pp.fold_runs_observed
        pp["super_additive"] = pp.pct_fold_runs_super_additive >= 50
        pp["feature_1"] = [k[0] for k in pp.key]; pp["feature_2"] = [k[1] for k in pp.key]
        pp[["feature_1", "feature_2", "interaction_information_mean", "fold_runs_observed", "pct_fold_runs_super_additive", "super_additive"]].to_csv(
            os.path.join(HERE, f"pair_interaction_{dt}.tsv"), sep="\t", index=False, float_format="%.3f")
        per_pair[dt] = {k: (a, b, bool(c)) for k, a, b, c in zip(pp.key, pp.interaction_information_mean, pp.pct_fold_runs_super_additive, pp.super_additive)}
        pp_tables = globals().setdefault("PP", {}); pp_tables[dt] = pp
    pd.DataFrame(summ).to_csv(os.path.join(HERE, "interaction_calibrated_summary.tsv"), sep="\t", index=False)

    cols = ["interaction_information_mean", "pct_fold_runs_super_additive", "super_additive"]
    def lookup(dt, a, b):
        k = tuple(sorted((a, b)))
        return list(per_pair[dt].get(k, (np.nan, np.nan, np.nan)))
    # stable pairs
    fin = pd.read_csv(os.path.join(REPO, "R2_classification_performance", "final_synthetic_features_table.csv"))
    vals = [lookup(r["Data Type"], r["Feature1"], r["Feature2"]) for _, r in fin.iterrows()]
    st = fin[["Data Type", "Feature1", "Feature2", "Selection %"]].copy(); st[cols] = vals
    st.to_csv(os.path.join(HERE, "stable_pairs_interaction.tsv"), sep="\t", index=False, float_format="%.3f")
    numbers.append({"section": "R5", "claim": "Stable pairs that are super-additive (>= 50 % of their fold-runs above the null threshold)",
                    "value": f"{int((st.super_additive == True).sum())}/{len(st)} (" + ", ".join(f"{k} {int((g.super_additive == True).sum())}/{len(g)}" for k, g in st.groupby('Data Type')) + ")",
                    "source": "R5 stable_pairs_interaction.tsv"})
    for a, b in (("Parvimonas_micra", "Peptostreptococcus_stomatis"), ("Fusobacterium_nucleatum", "Peptostreptococcus_stomatis"),
                 ("Dialister_pneumosintes", "Parvimonas_micra"), ("Fusobacterium_nucleatum", "Parvimonas_micra")):
        v = lookup("taxa", a, b)
        numbers.append({"section": "R5", "claim": f"{a} + {b} (taxa): mean interaction information; % fold-runs super-additive", "value": f"{v[0]:.0f} bits; {v[1]:.0f} %", "source": "R5 pair_interaction_taxa.tsv"})
    _pt = PP["taxa"]; top = _pt[_pt.fold_runs_observed >= 84].sort_values("interaction_information_mean", ascending=False).head(5).set_index("key")
    numbers.append({"section": "R5", "claim": "Taxa pairs with the highest mean interaction information (observed in >= 84 of 140 fold-runs)",
                    "value": "; ".join(f"{k[0]} + {k[1]} ({r.interaction_information_mean:.0f} bits, {r.pct_fold_runs_super_additive:.0f} %)" for k, r in top.iterrows()), "source": "R5 pair_interaction_taxa.tsv"})
    pd.DataFrame(numbers).to_csv(os.path.join(HERE, "r5_numbers_interaction.tsv"), sep="\t", index=False)
    # supplementary tables: add the columns
    p2 = os.path.join(MS, "SuppTable2_top30_synergistic_pairs.tsv"); s2 = pd.read_csv(p2, sep="\t").drop(columns=cols, errors="ignore")
    s2[cols] = [lookup(r.data_type, r.base_feature, r.contributing_feature) for r in s2.itertuples()]
    s2[cols[:2]] = s2[cols[:2]].astype(float).round(2); s2.to_csv(p2, sep="\t", index=False)
    p4 = os.path.join(MS, "SuppTable4_stable_pairs.tsv"); s4 = pd.read_csv(p4, sep="\t").drop(columns=cols, errors="ignore")
    s4[cols] = [lookup(r["Data Type"], r["Feature1"], r["Feature2"]) for _, r in s4.iterrows()]
    s4[cols[:2]] = s4[cols[:2]].astype(float).round(2); s4.to_csv(p4, sep="\t", index=False)
    print(pd.DataFrame(summ).round(3).to_string(index=False)); print(pd.DataFrame(numbers)[["claim", "value"]].to_string(index=False))


if __name__ == "__main__":
    main()
