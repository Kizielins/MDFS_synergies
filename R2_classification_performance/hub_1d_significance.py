#!/usr/bin/env python3
"""1D (univariate) MDFS significance of the hub species in every LOCO fold and MDFS run
(Results 2: "P. micra was individually significant in 1D MDFS in all 14 LOCO folds").

The 1D p-values of the LOCO runs are written by main_code.ipynb only since the revision
(results/<dt>/mdfs_1d_<dt>.csv). When that file is absent, the p-values are reconstructed from
the saved 1D IGs of the top-2 % pairs (mdfs_ig_info_<dt>.csv): under H0, 2*ln(2)*IG ~ chi2 with
df = response.divisions * divisions * (divisions+1)^(dimensions-1) = 1 for the default 1D MDFS
settings (MDFS::ComputePValue), and the per-feature p-value is Bonferroni-corrected over the number
of training features of the fold (more conservative than MDFS's Holm correction, which needs the
IGs of all features). Output: hub_1d_significance.tsv
Usage: python3 hub_1d_significance.py [feature ...]   (default: Parvimonas_micra Fusobacterium_nucleatum)
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RES = os.path.join(REPO, "M2_loco_crc_pipeline", "results")
PREP = os.path.join(REPO, "M1_data_preparation", "prepared_data")
FEATURES = sys.argv[1:] or ["Parvimonas_micra", "Fusobacterium_nucleatum"]
ALL_COHORT_STEMS = [
    'Public_study__FengQ_2015', 'Public_study__GuptaA_2019', 'Public_study__LiuNN_2022',
    'Public_study__ThomasAM_2018b', 'Public_study__VogtmannE_2016', 'Public_study__WirbelJ_2018',
    'Public_study__YachidaS_2019', 'Public_study__YangJ_2020', 'Public_study__YangY_2021',
    'Public_study__YuJ_2015', 'Public_study__ZellerG_2014', 'This_study_cohort1__ONCOBIOME_AtezoTRIBE',
    'This_study_cohort2__ONCOBIOME_COLOBIOME', 'This_study_cohort3__ONCOBIOME_IIGM_CZ',
    'This_study_cohort4__ONCOBIOME_IIGM_IT', 'This_study_cohort5__NSHII', 'This_study_cohort6__IIGM_TU'
]


def n_training_features(dt, test_cohort):
    cols = None
    for c in ALL_COHORT_STEMS:
        if c == test_cohort:
            continue
        f = os.path.join(PREP, f"X_{c}_{dt}.tsv")
        if os.path.exists(f):
            s = set(pd.read_csv(f, sep="\t", nrows=0, index_col=0).columns)
            cols = s if cols is None else cols & s
    return len(cols)


rows = []
for dt in ("taxa", "combined"):
    p1 = os.path.join(RES, dt, f"mdfs_1d_{dt}.csv")
    if os.path.exists(p1):                       # exact p-values saved by the pipeline
        d = pd.read_csv(p1)
        d = d[d["Feature"].isin(FEATURES)]
        for (f, fold, co), g in d.groupby(["Feature", "fold", "test_cohort"]):
            rows.append({"data_type": dt, "feature": f, "fold": fold, "test_cohort": co, "source": "mdfs_1d (MDFS Holm-adjusted p)",
                         "n_runs": len(g), "IG_min": g["IG"].min(), "IG_max": g["IG"].max(),
                         "p_adj_max": g["adjusted_p_value"].max(), "significant_all_runs": bool((g["adjusted_p_value"] < 0.05).all())})
        continue
    d = pd.read_csv(os.path.join(RES, dt, f"mdfs_ig_info_{dt}.csv"))
    nfeat = {co: n_training_features(dt, co) for co in d["test_cohort"].unique()}
    for f in FEATURES:
        a = d[d["base"] == f][["fold", "test_cohort", "run", "base_IG_1D"]].rename(columns={"base_IG_1D": "IG"})
        b = d[d["contributing"] == f][["fold", "test_cohort", "run", "contributing_IG_1D"]].rename(columns={"contributing_IG_1D": "IG"})
        ig = pd.concat([a, b]).drop_duplicates(["fold", "run"])
        for (fold, co), g in ig.groupby(["fold", "test_cohort"]):
            p = chi2.sf(2 * np.log(2) * g["IG"].values, df=1)
            p_adj = np.minimum(1.0, p * nfeat[co])
            rows.append({"data_type": dt, "feature": f, "fold": fold, "test_cohort": co,
                         "source": f"reconstructed from IG (chi2 df=1, Bonferroni x {nfeat[co]} features)",
                         "n_runs": len(g), "IG_min": g["IG"].min(), "IG_max": g["IG"].max(),
                         "p_adj_max": p_adj.max(), "significant_all_runs": bool((p_adj < 0.05).all())})
out = pd.DataFrame(rows).sort_values(["data_type", "feature", "fold"])
out.to_csv(os.path.join(HERE, "hub_1d_significance.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(out.to_string(index=False))
print("\nSummary (folds with the feature 1D-significant in every run it was observed / folds observed / 14):")
print(out.groupby(["data_type", "feature"]).agg(folds_sig=("significant_all_runs", "sum"), folds_obs=("fold", "nunique"),
                                                 runs_min=("n_runs", "min"), IG_min=("IG_min", "min")).to_string())
