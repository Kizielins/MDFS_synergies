"""Every number the manuscript edits need, read from the repository's final output files."""
import glob
import os
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

REPO = "/Users/kizie/Library/Mobile Documents/com~apple~CloudDocs/OneDrive-UJ/Sonata2/mdfs_combine_taxonomy_and_functions/CRC_segata_pooling/framework_repository"
MS = os.path.join(REPO, "manuscript")
M2 = os.environ.get("M2_RESULTS", os.path.join(REPO, "M2_loco_crc_pipeline", "results"))
R1, R2, R3, R4, R5 = (os.path.join(REPO, d) for d in ("R1_synergistic_pairs", "R2_classification_performance", "R3_pair_transformations", "R4_broader_application", "R5_reviewer_checks"))
SP = {"Fusobacterium_nucleatum": "F. nucleatum", "Parvimonas_micra": "P. micra", "Peptostreptococcus_stomatis": "P. stomatis",
      "Dialister_pneumosintes": "D. pneumosintes", "Anaerococcus_obesiensis": "A. obesiensis", "Gemella_morbillorum": "G. morbillorum",
      "Roseburia_intestinalis": "R. intestinalis", "Eikenella_corrodens": "E. corrodens", "Blautia_hydrogenotrophica": "B. hydrogenotrophica"}


def sp(name):
    if name in SP:
        return SP[name]
    if ":" in name:                      # pathway: "P163-PWY:_L-lysine_fermentation..." -> "L-lysine fermentation to acetate and butanoate (P163-PWY)"
        pid, desc = name.split(":", 1)
        return f"{desc.strip('_').replace('_', ' ')} ({pid})"
    parts = name.split("_")
    if len(parts) >= 2 and parts[0][0].isupper() and parts[1][0].islower():
        return f"{parts[0][0]}. {' '.join(parts[1:])}"
    return name.replace("_", " ")


def fmt(x, nd=0):
    return f"{x:,.{nd}f}"


def load():
    N = {}
    # ---------------- R2: AUCs ----------------
    perf = pd.read_csv(os.path.join(R2, "loco_performance_summary.tsv"), sep="\t")
    key = {"Feature Only": "raw", "Synthetic Only": "syn", "Combined": "both"}
    for _, r in perf.iterrows():
        N[f"auc_{r['Mode'].lower()}_{key[r['Feature Type']]}"] = r["Mean AUC"]
        N[f"wauc_{r['Mode'].lower()}_{key[r['Feature Type']]}"] = r["Weighted Mean AUC"]
    for m in ("taxa", "functions", "combined"):
        v = [N[f"auc_{m}_{k}"] for k in ("raw", "syn", "both")]
        N[f"auc_{m}_range"] = f"{min(v):.2f}" if f"{min(v):.2f}" == f"{max(v):.2f}" else f"{min(v):.2f}-{max(v):.2f}"
    # ---------------- R2: stable pairs ----------------
    fin = pd.read_csv(os.path.join(R2, "final_synthetic_features_table.csv"))
    N["n_stable"] = len(fin)
    N["n_stable_by_dt"] = fin.groupby("Data Type").size().to_dict()
    N["n_stable_zero"] = int(((fin["Feature1_Selection%"] == 0) | (fin["Feature2_Selection%"] == 0)).sum())
    base = pd.read_csv(os.path.join(M2, "combined", "base_features_combined.csv")).set_index("Feature")["Selection %"]
    for f, k in (("Parvimonas_micra", "pm"), ("Fusobacterium_nucleatum", "fn")):
        N[f"{k}_raw_freq"] = float(base.get(f, 0.0))
        sub = fin[(fin["Data Type"] == "combined") & ((fin.Feature1 == f) | (fin.Feature2 == f))]
        as_base = sum(((r.Feature1 == f) and (r["Feature1_Selection%"] <= r["Feature2_Selection%"])) or ((r.Feature2 == f) and (r["Feature2_Selection%"] < r["Feature1_Selection%"])) for _, r in sub.iterrows())
        N[f"{k}_stable_combined"] = len(sub); N[f"{k}_as_base"] = as_base; N[f"{k}_as_contrib"] = len(sub) - as_base
        N[f"{k}_stable_all"] = int(((fin.Feature1 == f) | (fin.Feature2 == f)).sum())
    meta = pd.read_csv(os.path.join(R2, "meta_results_all_pairs.csv"))
    N["n_meta_sig"] = int((meta.q_meta_fdr < 0.05).sum()); N["n_meta"] = len(meta)
    cons = meta.direction.apply(lambda d: max(int(d.split("/")[0].lstrip("+")), int(d.split("/")[1].lstrip("-"))))
    N["n_cohorts_meta"] = int(meta.n_cohorts.max())
    N["dir14"] = int((cons >= 14).sum()); N["dir15"] = int((cons == 15).sum()); N["dir13"] = int((cons == 13).sum()); N["dir12"] = int((cons == 12).sum()); N["dir_min"] = int(cons.min())
    ct = pd.read_csv(os.path.join(R2, "cohort_level_tests.csv"))
    N["cohort_tests_sig"] = int((ct.p_value_fdr < 0.05).sum()); N["cohort_tests"] = len(ct); N["cohort_tests_pct"] = 100 * (ct.p_value_fdr < 0.05).mean()
    # ---------------- R3 ----------------
    s3 = pd.read_csv(os.path.join(MS, "SuppTable3_transformation_preference.tsv"), sep="\t"); s3["types"] = s3.Types_Present.astype(str)
    m = fin.merge(s3, on=["Data Type", "Feature1", "Feature2"], how="left")
    m2 = fin.merge(s3.rename(columns={"Feature1": "Feature2", "Feature2": "Feature1"}), on=["Data Type", "Feature1", "Feature2"], how="left")
    m["types"] = m["types"].fillna(m2["types"]); m["Num_Types"] = m["Num_Types"].fillna(m2["Num_Types"])
    N["lr_only"] = int((m.types == "['LR']").sum()); N["gm_only"] = int((m.types == "['GM']").sum()); N["both"] = int((m.Num_Types == 2).sum())
    N["gm_based"] = N["gm_only"] + N["both"]; N["gm_based_pct"] = 100 * N["gm_based"] / N["n_stable"]
    N["n_supp3"] = len(s3)
    cand = pd.read_csv(os.path.join(R3, "figure4_candidates.tsv"), sep="\t")
    # the nine Figure 4 pairs in panel order are the caption printed by the notebook; recover from the notebook output
    import json
    nb = json.load(open(os.path.join(R3, "transformation_preference.ipynb")))
    fig4 = []
    for c in nb["cells"]:
        for o in c.get("outputs", []):
            t = "".join(o.get("text", [])) if "text" in o else ""
            if "Selected for Figure 4:" in t:
                block = t.split("Selected for Figure 4:")[1].strip().splitlines()[1:10]
                for line in block:
                    parts = line.split()
                    fig4.append({"data_type": parts[0], "pair": parts[1], "sel": float(parts[2]), "types": parts[3], "lr": float(parts[4]), "gm": float(parts[5])})
    N["fig4"] = fig4
    # ---------------- R1 ----------------
    N["wilcoxon"] = {}
    N["consistency"] = {}
    N["fn_partners"] = {}
    for dt in ("taxa", "functions", "combined"):
        g = pd.read_csv(os.path.join(R1, f"ig_gain_summary_{dt}.csv"))
        N["wilcoxon"][dt] = wilcoxon(g["total_IG_2D_mean"], g["base_IG_1D_mean"], alternative="greater").pvalue
        g["n_obs"] = g.pair_ig_count.str.split("/").str[0].astype(int)
        g["key"] = [tuple(sorted((a, b))) for a, b in zip(g.base_feature, g.contributing_feature)]
        best = g.groupby("key").n_obs.max().sort_values(ascending=False)
        N["consistency"][dt] = [(k, int(v)) for k, v in best.head(5).items()]
        # top 10 unordered pairs by mean total 2D IG (the Supplementary Table 1 ranking)
        top = g.sort_values("total_IG_2D_mean", ascending=False).drop_duplicates("key").head(10)
        N["fn_partners"][dt] = []
        for k in top.key:
            if "Fusobacterium_nucleatum" in k:
                N["fn_partners"][dt].append(k[0] if k[1] == "Fusobacterium_nucleatum" else k[1])
        N[f"fn_in_top10_{dt}"] = len(N["fn_partners"][dt])
        N.setdefault("top10", {})[dt] = list(top.key)
    s5 = pd.read_csv(os.path.join(MS, "SuppTable5_seed_stability.tsv"), sep="\t")
    N["jac"] = {dt: (float(s5[(s5.data_type == dt) & (s5.N == 10) & (s5.k == 6)].jaccard_mean.iloc[0]), float(s5[(s5.data_type == dt) & (s5.N == 2) & (s5.k == 2)].jaccard_mean.iloc[0])) for dt in ("taxa", "functions", "combined")}
    fb = pd.read_csv(os.path.join(R1, "figure1b_all_pairs", "figure1b_all_directed_summary.tsv"), sep="\t")
    N["fig1b"] = fb.set_index("data_type").to_dict("index")
    hub = pd.read_csv(os.path.join(R2, "hub_1d_significance.tsv"), sep="\t")
    for f, k in (("Parvimonas_micra", "pm"), ("Fusobacterium_nucleatum", "fn")):
        for dt in ("taxa", "combined"):
            h = hub[(hub.feature == f) & (hub.data_type == dt)]
            N[f"{k}_1d_folds_{dt}"] = int(h.significant_all_runs.sum()); N[f"{k}_1d_source"] = h.source.iloc[0] if len(h) else ""
    # ---------------- R4 ----------------
    t1 = pd.read_csv(os.path.join(MS, "Table1_validation_cohorts.tsv"), sep="\t")
    an = t1[t1.status == "analysed"]
    N["n_meta_cohorts"] = len(t1); N["n_excluded_crc"] = int(t1.status.str.startswith("excluded").sum()); N["n_noncrc"] = len(t1) - N["n_excluded_crc"]
    N["n_analysed"] = len(an); N["n_with_pairs"] = int((an.n_synergistic_pairs > 0).sum())
    N["cohorts_with_pairs"] = an[an.n_synergistic_pairs > 0].sort_values("n_synergistic_pairs", ascending=False)[["cohort", "n_synergistic_pairs", "n_pairs_gain_ge_50pct", "conditions"]].to_dict("records")
    N["cohorts_zero"] = list(an[an.n_synergistic_pairs == 0].cohort)
    N["table1"] = t1.set_index("cohort").to_dict("index")
    syn = pd.read_csv(os.path.join(R4, "synergy_by_cohort.tsv"), sep="\t")
    syn["pair_key"] = [" + ".join(sorted((str(a), str(b)))) for a, b in zip(syn.base, syn.contributing)]
    synu = syn.sort_values("total_IG_mean", ascending=False).drop_duplicates(["cohort", "pair_key"])
    N["n_pairs_total"] = int(synu.pair_key.nunique())
    N["n_pairs_100"] = int(synu.loc[synu.ig_gain_pct >= 100, "pair_key"].nunique()); N["pct_pairs_100"] = 100 * N["n_pairs_100"] / N["n_pairs_total"]
    N["share100"] = {co: 100 * (synu[synu.cohort == co].ig_gain_pct >= 100).mean() for co in synu.cohort.unique()}
    N["n_pairs_cohort"] = synu.groupby("cohort").pair_key.nunique().to_dict()
    N["n_pairs50_cohort"] = synu[synu.ig_gain_pct >= 50].groupby("cohort").pair_key.nunique().to_dict()
    emg = synu[(synu.base_IG_1D_mean == 0) | (synu.contributing_IG_1D_mean == 0)]
    N["n_emergent"] = int(emg.pair_key.nunique()); N["emergent_cohorts"] = sorted(emg.cohort.unique())
    def pair(co, a, b):
        d = syn[(syn.cohort == co) & (((syn.base == a) & (syn.contributing == b)) | ((syn.base == b) & (syn.contributing == a)))]
        d = d[d.base == a]
        return None if not len(d) else dict(base1d=float(d.base_IG_1D_mean.iloc[0]), added=float(d.IG_2D_added_mean.iloc[0]), total=float(d.total_IG_mean.iloc[0]), runs=int(d.n_runs.iloc[0]))
    N["eramulus"] = pair("HallAB_2017", "Eubacterium_ramulus", "Ruminococcus_torques")
    N["camalon"] = pair("QinN_2014", "Citrobacter_amalonaticus", "Streptococcus_anginosus_group")
    mc = synu[synu.cohort == "MetaCardis_2020_a"].sort_values("total_IG_mean", ascending=False).head(2)
    N["metacardis_top"] = [(r.pair_key, float(r.total_IG_mean)) for r in mc.itertuples()]
    gem = synu[synu.pair_key.str.contains("Gemmiger_formicilis") & synu.pair_key.str.contains("PWY-6609")]
    N["gemmiger_cohorts"] = int(gem.cohort.nunique())
    nco = synu.groupby("pair_key").cohort.nunique()
    N["max_cohorts"] = int(nco.max()); N["n_pairs_ge2"] = int((nco >= 2).sum()); N["n_pairs_maxco"] = int((nco == nco.max()).sum())
    top3 = synu[synu.pair_key.isin(nco[nco == nco.max()].index)]
    N["maxco_cohorts"] = sorted(top3.cohort.unique())
    N["maxco_strep_share"] = 100 * np.mean([("Streptococcus" in k) for k in nco[nco == nco.max()].index]) if N["n_pairs_maxco"] else 0
    s4 = pd.read_csv(os.path.join(MS, "SuppTable4_conserved_pairs.tsv"), sep="\t")
    N["n_conserved"] = len(s4); N["conserved_diseases"] = s4.groupby("disease").size().to_dict()
    N["conserved_taxa_taxa"] = int(sum(":" not in k for k in s4.pair_key)); N["conserved_taxa_fun"] = len(s4) - N["conserved_taxa_taxa"]
    N["conserved_same_orient"] = int(s4.same_orientation_in_all_cohorts.sum()) if "same_orientation_in_all_cohorts" in s4 else None
    N["conserved_top"] = s4.sort_values("mean_total_ig", ascending=False).head(3)[["pair_key", "mean_total_ig", "mean_ig_gain_pct"]].to_dict("records")
    N["conserved_common_feature"] = None
    from collections import Counter
    cnt = Counter(f for k in s4.pair_key for f in k.split(" + "))
    if cnt and cnt.most_common(1)[0][1] == len(s4):
        N["conserved_common_feature"] = cnt.most_common(1)[0][0]
    # ---------------- R5 (optional) ----------------
    N["r5"] = {}
    for f in glob.glob(os.path.join(R5, "r5_numbers*.tsv")):
        for _, r in pd.read_csv(f, sep="\t").iterrows():
            N["r5"][r["claim"]] = r["value"]
    for f in ("interaction_information.tsv", "stable_pairs_lrt.tsv", "collinearity.tsv", "negative_controls.tsv", "discretisation_sensitivity.tsv", "runtime.tsv", "sensitivity_auc.tsv"):
        p = os.path.join(R5, f)
        N["r5_" + f.split(".")[0]] = pd.read_csv(p, sep="\t") if os.path.exists(p) else None
    return N


if __name__ == "__main__":
    N = load()
    for k, v in N.items():
        if not k.startswith("r5_") and k not in ("table1",):
            print(k, "=", v if not isinstance(v, float) else round(v, 3))
