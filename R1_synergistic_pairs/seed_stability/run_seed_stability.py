#!/usr/bin/env python3
"""
MDFS seed-stability analysis for the LOCO discovery arm.

For every LOCO fold (same training sets as main_code.ipynb: all other cohorts,
sorted common features) and data type, run run_mdfs.R with N_SEEDS different
seeds and store the top 5% unordered pairs of each run, scored by the larger
total_IG of the two directions (the top 2% used by main_code is a prefix of this).

Output: pairs/<data_type>_fold<NN>_seed<SEED>.tsv.gz  (rank, feature_1, feature_2, score, n_unique_pairs)
No Random Forest is trained; nothing in framework_repository is modified.
"""
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PREP = os.path.join(REPO, "M1_data_preparation", "prepared_data")
RSCRIPT = os.path.join(REPO, "run_mdfs.R")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pairs")
DATA_TYPES = ["taxa", "functions", "combined"]
SEEDS = list(range(42, 62))      # 20 seeds; 42..51 = the seeds used by main_code.ipynb
KEEP_FRACTION = 0.05
N_WORKERS = 12

ALL_COHORT_STEMS = [
    'Public_study__FengQ_2015', 'Public_study__GuptaA_2019', 'Public_study__LiuNN_2022',
    'Public_study__ThomasAM_2018b', 'Public_study__VogtmannE_2016', 'Public_study__WirbelJ_2018',
    'Public_study__YachidaS_2019', 'Public_study__YangJ_2020', 'Public_study__YangY_2021',
    'Public_study__YuJ_2015', 'Public_study__ZellerG_2014', 'This_study_cohort1__ONCOBIOME_AtezoTRIBE',
    'This_study_cohort2__ONCOBIOME_COLOBIOME', 'This_study_cohort3__ONCOBIOME_IIGM_CZ',
    'This_study_cohort4__ONCOBIOME_IIGM_IT', 'This_study_cohort5__NSHII', 'This_study_cohort6__IIGM_TU'
]


def training_set(data_type, test_cohort):
    """Replicates main_code.ipynb: all other cohorts, sorted common features, concatenated."""
    Xs, ys = [], []
    for c in ALL_COHORT_STEMS:
        if c == test_cohort:
            continue
        xf = os.path.join(PREP, f"X_{c}_{data_type}.tsv")
        yf = os.path.join(PREP, f"y_{c}_{data_type}.tsv")
        if not (os.path.exists(xf) and os.path.exists(yf)):
            continue
        X = pd.read_csv(xf, sep="\t", index_col=0)
        y = pd.read_csv(yf, sep="\t", index_col=0).squeeze().astype(int)
        common = X.index.intersection(y.index)
        Xs.append(X.loc[common])
        ys.append(y.loc[common])
    feats = sorted(set.intersection(*[set(x.columns) for x in Xs]))
    return pd.concat([x[feats] for x in Xs]), pd.concat(ys)


def run_fold(data_type, fold, test_cohort):
    X, y = training_set(data_type, test_cohort)
    done = 0
    with tempfile.TemporaryDirectory(prefix="seedstab_") as tmp:
        xp, yp = os.path.join(tmp, "X.tsv"), os.path.join(tmp, "y.tsv")
        X.to_csv(xp, sep="\t")
        y.rename("label").to_csv(yp, sep="\t", header=True)
        for seed in SEEDS:
            out = os.path.join(OUT, f"{data_type}_fold{fold:02d}_seed{seed}.tsv.gz")
            if os.path.exists(out):
                done += 1
                continue
            p2d, pu, p1d = (os.path.join(tmp, n) for n in ("2d.tsv", "u.txt", "1d.tsv"))
            r = subprocess.run(["Rscript", RSCRIPT, xp, yp, p2d, pu, p1d, str(seed)],
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"{data_type} fold {fold} seed {seed}: {r.stderr[-500:]}")
            pairs = pd.read_csv(p2d, sep="\t")
            score = {}
            for b, c, t in zip(pairs["base"].astype(str).str.replace(" ", "_"),
                               pairs["contributing"].astype(str).str.replace(" ", "_"),
                               pairs["total_IG"].astype(float)):
                k = tuple(sorted((b, c)))
                score[k] = max(score.get(k, 0.0), t)
            ranked = sorted(score.items(), key=lambda kv: -kv[1])
            n_keep = max(1, int(len(ranked) * KEEP_FRACTION))
            df = pd.DataFrame([(i + 1, a, b, s) for i, ((a, b), s) in enumerate(ranked[:n_keep])],
                              columns=["rank", "feature_1", "feature_2", "score"])
            df["n_unique_pairs"] = len(ranked)
            df.to_csv(out + ".tmp", sep="\t", index=False, compression="gzip")
            os.replace(out + ".tmp", out)
            done += 1
    return data_type, fold, test_cohort, X.shape, done


def main():
    os.makedirs(OUT, exist_ok=True)
    tasks = []
    for dt in DATA_TYPES:
        # Same folds / test cohorts as the saved LOCO results
        folds = pd.read_csv(os.path.join(REPO, f"results_{dt}", f"mdfs_ig_info_{dt}.csv"),
                            usecols=["fold", "test_cohort"]).drop_duplicates().sort_values("fold")
        tasks += [(dt, int(f), c) for f, c in zip(folds["fold"], folds["test_cohort"])]
    print(f"{len(tasks)} fold tasks x {len(SEEDS)} seeds", flush=True)
    with ProcessPoolExecutor(N_WORKERS) as ex:
        futs = [ex.submit(run_fold, *t) for t in tasks]
        for i, f in enumerate(as_completed(futs), 1):
            dt, fold, tc, shape, done = f.result()
            print(f"[{i}/{len(tasks)}] {dt} fold {fold} ({tc}) train {shape}: {done} seeds", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
