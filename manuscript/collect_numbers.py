#!/usr/bin/env python3
"""Collect every number quoted in the Results / Methods text from the final output files and write
RESULTS_NUMBERS.md (one table per manuscript section, each row = claim, value, source file).
Run after run_all.sh:  python3 manuscript/collect_numbers.py
"""
import glob
import os
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
M2 = os.path.join(REPO, "M2_loco_crc_pipeline", "results")
R1, R2, R3, R4 = (os.path.join(REPO, d) for d in ("R1_synergistic_pairs", "R2_classification_performance",
                                                   "R3_pair_transformations", "R4_broader_application"))
rows = []
def add(section, claim, value, source):
    rows.append((section, claim, value, source))

# ---------------- Results 1: synergistic pairs ----------------
for dt in ("taxa", "functions", "combined"):
    g = pd.read_csv(os.path.join(R1, f"ig_gain_summary_{dt}.csv"))
    p = wilcoxon(g["total_IG_2D_mean"], g["base_IG_1D_mean"], alternative="greater").pvalue
    add("R1", f"Wilcoxon 2D vs 1D IG, {dt} (directed pairs in the per-run top 2 %)", f"n = {len(g):,}, p = {p:.2g}",
        f"R1_synergistic_pairs/ig_gain_summary_{dt}.csv")
    g["n_obs"] = g["pair_ig_count"].str.split("/").str[0].astype(int)
    full = g[g["n_obs"] == 140]
    add("R1", f"Directed pairs in the top 2 % of all 140 fold-runs, {dt}", f"{len(full)}: " + "; ".join(
        f"{a} -> {b}" for a, b in zip(full["base_feature"], full["contributing_feature"])), f"ig_gain_summary_{dt}.csv")
    top = g.sort_values("ig_difference_raw", ascending=False).head(10)
    n_fn = sum(("Fusobacterium_nucleatum" in a) or ("Fusobacterium_nucleatum" in b) for a, b in zip(top["base_feature"], top["contributing_feature"]))
    add("R1", f"Top-10 pairs by raw 2D gain containing F. nucleatum, {dt}", f"{n_fn}/10", f"ig_gain_summary_{dt}.csv")
    for a, b in (("Fusobacterium_nucleatum", "Peptostreptococcus_stomatis"), ("Fusobacterium_nucleatum", "Parvimonas_micra"),
                 ("Fusobacterium_nucleatum", "Dialister_pneumosintes"), ("Fusobacterium_nucleatum", "Anaerococcus_obesiensis"),
                 ("Parvimonas_micra", "Dialister_pneumosintes"), ("Parvimonas_micra", "Peptostreptococcus_stomatis")):
        sel = g[((g.base_feature == a) & (g.contributing_feature == b)) | ((g.base_feature == b) & (g.contributing_feature == a))]
        if len(sel):
            add("R1", f"{a} + {b}, {dt}: fold-runs in top 2 % (max over the two directions)", f"{sel['n_obs'].max()}/140", f"ig_gain_summary_{dt}.csv")
s5 = pd.read_csv(os.path.join(HERE, "SuppTable1_seed_stability.tsv"), sep="\t")
for dt in ("taxa", "functions", "combined"):
    r = s5[(s5.data_type == dt) & (s5.N == 10) & (s5.k == 6)].iloc[0]
    r2 = s5[(s5.data_type == dt) & (s5.N == 2) & (s5.k == 2)].iloc[0]
    add("R1/Methods", f"Seed stability {dt}: Jaccard of consensus sets from two disjoint seed blocks",
        f"10 runs, >= 6: {r.jaccard_mean:.2f} (min {r.jaccard_min:.2f}); 2 runs, both: {r2.jaccard_mean:.2f}", "SuppTable1_seed_stability.tsv")

# ---------------- Results 2: classification ----------------
perf = pd.read_csv(os.path.join(R2, "loco_performance_summary.tsv"), sep="\t")
for _, r in perf.iterrows():
    add("R2", f"Mean AUC over 14 test cohorts, {r['Mode']} / {r['Feature Type']}", f"{r['Mean AUC']:.3f} (weighted {r['Weighted Mean AUC']:.3f})",
        "R2_classification_performance/loco_performance_summary.tsv (from M2 results loco_results_<dt>.csv)")
for dt in ("taxa", "functions", "combined"):
    lr = pd.read_csv(os.path.join(M2, dt, f"loco_results_{dt}.csv"))
    add("R2", f"LOCO test cohorts (folds), {dt}", f"{len(lr) - 2}", f"M2_loco_crc_pipeline/results/{dt}/loco_results_{dt}.csv")
fin = pd.read_csv(os.path.join(R2, "final_synthetic_features_table.csv"))
add("R2", "Stable pairs (>= 80 % of folds, frequency above both constituents)", f"{len(fin)} (" + ", ".join(f"{k} {v}" for k, v in fin.groupby("Data Type").size().items()) + ")",
    "R2_classification_performance/final_synthetic_features_table.csv")
zero = ((fin["Feature1_Selection%"] == 0) | (fin["Feature2_Selection%"] == 0)).sum()
add("R2", "Stable pairs with >= 1 constituent at 0 % individual frequency", f"{zero}/{len(fin)}", "final_synthetic_features_table.csv")
base = pd.read_csv(os.path.join(M2, "combined", "base_features_combined.csv")).set_index("Feature")["Selection %"]
for f in ("Parvimonas_micra", "Fusobacterium_nucleatum"):
    add("R2", f"{f}: individual selection frequency in the combined model", f"{base.get(f, 0.0):.0f} %", "M2 results/combined/base_features_combined.csv")
    sub = fin[(fin["Data Type"] == "combined") & ((fin.Feature1 == f) | (fin.Feature2 == f))]
    # base = the constituent with the lower individual frequency (as in Figure 3)
    as_base = sum(((r.Feature1 == f) and (r["Feature1_Selection%"] <= r["Feature2_Selection%"])) or ((r.Feature2 == f) and (r["Feature2_Selection%"] < r["Feature1_Selection%"])) for _, r in sub.iterrows())
    add("R2", f"{f}: stable combined-model pairs (as base = lower-frequency constituent / as contributor)", f"{len(sub)} ({as_base} base / {len(sub) - as_base} contributor)", "final_synthetic_features_table.csv")
meta = pd.read_csv(os.path.join(R2, "meta_results_all_pairs.csv"))
add("R2", "Pairs meta-significant (Stouffer, BH q < 0.05)", f"{(meta.q_meta_fdr < 0.05).sum()}/{len(meta)}", "R2_classification_performance/meta_results_all_pairs.csv")
add("R2", "Cohorts per pair in the meta-analysis", ", ".join(map(str, sorted(meta.n_cohorts.unique()))), "meta_results_all_pairs.csv")
cons = meta.direction.apply(lambda d: max(int(d.split("/")[0].lstrip("+")), int(d.split("/")[1].lstrip("-"))))
add("R2", "Direction consistent in >= 14 of 15 cohorts (15/15; 13/15; 12/15)",
    f"{(cons >= 14).sum()} ({(cons == 15).sum()}; {(cons == 13).sum()}; {(cons == 12).sum()})", "meta_results_all_pairs.csv")
ct = os.path.join(R2, "cohort_level_tests.csv")
if os.path.exists(ct):
    c = pd.read_csv(ct)
    add("R2", "Cohort-pair tests significant (BH within pair, q < 0.05)", f"{(c.p_value_fdr < 0.05).sum()}/{len(c)} ({100 * (c.p_value_fdr < 0.05).mean():.1f} %)", "R2_classification_performance/cohort_level_tests.csv")

# ---------------- Results 3: transformations ----------------
s3 = pd.read_csv(os.path.join(HERE, "SuppTable6_transformation_preference.tsv"), sep="\t")
s3["types"] = s3["Types_Present"].astype(str)
LR_ONLY, GM_ONLY = "['LR']", "['GM']"
add("R3", "Pairs selected in >= 50 % of folds (Supp Table 6): LR only / GM only / both",
    f"{len(s3)}: {(s3.types == LR_ONLY).sum()} / {(s3.types == GM_ONLY).sum()} / {(s3.Num_Types == 2).sum()}", "SuppTable6_transformation_preference.tsv")
m = fin.merge(s3, on=["Data Type", "Feature1", "Feature2"], how="left")
m2 = fin.merge(s3.rename(columns={"Feature1": "Feature2", "Feature2": "Feature1"}), on=["Data Type", "Feature1", "Feature2"], how="left")
m["types"] = m["types"].fillna(m2["types"]); m["Num_Types"] = m["Num_Types"].fillna(m2["Num_Types"])
lr, gm, both = (m.types == LR_ONLY).sum(), (m.types == GM_ONLY).sum(), (m.Num_Types == 2).sum()
add("R3", f"Stable pairs ({len(fin)}): LR only / GM only / both (>= 50 % per transformation); GM-based = GM only + both",
    f"{lr} / {gm} / {both}; GM-based {gm + both}/{len(fin)} ({100 * (gm + both) / len(fin):.0f} %)", "final_synthetic_features_table.csv x SuppTable6")
for a, b in (("Blautia_wexlerae", "Parvimonas_micra"), ("Eikenella_corrodens", "Parvimonas_micra")):
    r = m[(m.Feature1 == a) & (m.Feature2 == b) | (m.Feature1 == b) & (m.Feature2 == a)]
    add("R3", f"{a} + {b}: stable? transformations at >= 50 %", f"{'yes' if len(r) else 'NOT stable'}; {r.types.iloc[0] if len(r) else s3[(s3.Feature1 == a) & (s3.Feature2 == b)].types.tolist()}", "final table / SuppTable6")

# ---------------- Results 4: broader application ----------------
t1 = pd.read_csv(os.path.join(HERE, "Table1_validation_cohorts.tsv"), sep="\t")
an = t1[t1.status == "analysed"]
add("R4", "Cohorts in metadata / excluded CRC / analysed / with >= 1 synergistic pair", f"{len(t1)} / {(t1.status.str.startswith('excluded')).sum()} / {len(an)} / {(an.n_synergistic_pairs > 0).sum()}", "Table1_validation_cohorts.tsv")
add("R4", "Analysed cohorts with 0 pairs", ", ".join(an[an.n_synergistic_pairs == 0].cohort), "Table1_validation_cohorts.tsv")
syn = pd.read_csv(os.path.join(R4, "synergy_by_cohort.tsv"), sep="\t")
syn["pair_key"] = [" + ".join(sorted((str(a), str(b)))) for a, b in zip(syn.base, syn.contributing)]
syn = syn.sort_values("total_IG_mean", ascending=False).drop_duplicates(["cohort", "pair_key"])
add("R4", "Unique (unordered) synergistic pairs over all cohorts", f"{syn.pair_key.nunique():,}", "R4_broader_application/synergy_by_cohort.tsv (unordered)")
add("R4", "Unique pairs with >= 100 % gain in >= 1 cohort", f"{syn.loc[syn.ig_gain_pct >= 100, 'pair_key'].nunique():,} ({100 * syn.loc[syn.ig_gain_pct >= 100, 'pair_key'].nunique() / syn.pair_key.nunique():.1f} %)", "synergy_by_cohort.tsv (unordered)")
for co in ("QinN_2014", "JieZ_2017", "MetaCardis_2020_a", "HallAB_2017"):
    s = syn[syn.cohort == co]
    add("R4", f"{co}: pairs; share with >= 100 % gain (all pairs); share among pairs with total IG >= 5 (Figure 5c filter)",
        f"{len(s):,}; {100 * (s.ig_gain_pct >= 100).mean():.1f} %; {100 * (s[s.total_IG_mean >= 5].ig_gain_pct >= 100).mean():.1f} %", "synergy_by_cohort.tsv (unordered)")
add("R4", "Pairs per cohort with >= 50 % gain (Figure 5b)", "; ".join(f"{k} {v:,}" for k, v in syn[syn.ig_gain_pct >= 50].groupby("cohort").size().sort_values(ascending=False).items()), "synergy_by_cohort.tsv (unordered)")
emg = syn[(syn.base_IG_1D_mean == 0) | (syn.contributing_IG_1D_mean == 0)]
add("R4", "Emergent pairs (one constituent with 1D IG = 0) and their cohorts", f"{emg.pair_key.nunique():,}; {', '.join(sorted(emg.cohort.unique()))}", "synergy_by_cohort.tsv (unordered)")
def pair(co, a, b):
    s = syn[(syn.cohort == co) & (syn.pair_key == " + ".join(sorted((a, b))))]
    if not len(s): return "absent"
    r = s.iloc[0]; d = pd.read_csv(os.path.join(R4, "synergy_by_cohort.tsv"), sep="\t")
    d = d[(d.cohort == co) & (((d.base == a) & (d.contributing == b)) | ((d.base == b) & (d.contributing == a)))]
    return "; ".join(f"{x.base} -> {x.contributing}: 1D {x.base_IG_1D_mean:.2f} + added {x.IG_2D_added_mean:.1f} = {x.total_IG_mean:.1f} ({x.n_runs}/10 runs, x{x.total_IG_mean / x.base_IG_1D_mean if x.base_IG_1D_mean else float('inf'):.0f})" for x in d.itertuples())
add("R4", "E. ramulus + R. torques, HallAB_2017", pair("HallAB_2017", "Eubacterium_ramulus", "Ruminococcus_torques"), "synergy_by_cohort.tsv")
add("R4", "C. amalonaticus + S. anginosus group, QinN_2014", pair("QinN_2014", "Citrobacter_amalonaticus", "Streptococcus_anginosus_group"), "synergy_by_cohort.tsv")
mc = syn[syn.cohort == "MetaCardis_2020_a"].sort_values("total_IG_mean", ascending=False).head(3)
add("R4", "MetaCardis_2020_a: top pairs by total IG", "; ".join(f"{r.pair_key} ({r.total_IG_mean:.1f})" for r in mc.itertuples()), "synergy_by_cohort.tsv")
gem = syn[syn.pair_key.str.contains("Gemmiger_formicilis") & syn.pair_key.str.contains("PWY-6609")]
add("R4", "Gemmiger formicilis + adenine/adenosine salvage III (PWY-6609): cohorts (>= 6/10 runs)", f"{gem.cohort.nunique()} ({', '.join(sorted(gem.cohort.unique()))})", "synergy_by_cohort.tsv")
nco = syn.groupby("pair_key").cohort.nunique()
add("R4", "Max number of cohorts sharing a pair; pairs in >= 2 / >= 3 cohorts", f"{nco.max()}; {(nco >= 2).sum():,} / {(nco >= 3).sum():,}", "synergy_by_cohort.tsv (unordered)")
s4 = pd.read_csv(os.path.join(HERE, "SuppTable7_conserved_pairs.tsv"), sep="\t")
add("R4", "Conserved pairs, >= 2 cohorts with the same disease label (Supp Table 7)", f"{len(s4)} (" + ", ".join(f"{k} {v}" for k, v in s4.groupby('disease').size().items()) + ")", "SuppTable7_conserved_pairs.tsv")

# ---------------- R5: reviewer checks (if run) ----------------
for f in sorted(glob.glob(os.path.join(REPO, "R5_reviewer_checks", "r5_numbers*.tsv"))):
    for _, r in pd.read_csv(f, sep="\t").iterrows():
        add("R5 reviewer checks", r["claim"], r["value"], r["source"])

out = ["# Numbers quoted in the manuscript and where they come from", "",
       f"Generated by `manuscript/collect_numbers.py` from the final output files (run {pd.Timestamp.now():%Y-%m-%d %H:%M}).",
       "Sections: R1-R4 = Results subsections in manuscript order. Compare each row with the text before submission.", ""]
for sec in dict.fromkeys(s for s, *_ in rows):
    out += [f"## {sec}", "", "| Claim | Value | Source |", "|---|---|---|"]
    out += [f"| {c} | {v} | `{s}` |" for sc, c, v, s in rows if sc == sec]
    out.append("")
open(os.path.join(HERE, "RESULTS_NUMBERS.md"), "w").write("\n".join(out))
print("\n".join(out))
