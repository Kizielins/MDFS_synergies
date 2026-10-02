#!/usr/bin/env python3
"""
Figure 1b input (all points): every directed pair (base, contributing) per LOCO fold and MDFS run
(same training sets and seeds 42..51 as main_code.ipynb), aggregated per fold. Checks that each run's
synergistic pairs (top 2%) equal the pairs saved by the LOCO pipeline (M2_loco_crc_pipeline/results).
Runtime: several hours (MDFS on every fold x 10 seeds, 6 processes).
Output: directed/<dt>_fold<NN>.tsv.gz, check_all_directed.tsv
Usage: python run_all_directed.py [dt:fold ...]   (arguments = test subset)
"""
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RESULTS = os.path.join(REPO, "M2_loco_crc_pipeline", "results")
PREP = os.path.join(REPO, "M1_data_preparation", "prepared_data")
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

SEEDS = list(range(42, 52))
OUT = os.path.join(HERE, "directed")


def run_fold(dt, fold, test_cohort):
    X, y = training_set(dt, test_cohort)
    saved = pd.read_csv(os.path.join(RESULTS, dt, f"mdfs_ig_info_{dt}.csv"),
                        usecols=["fold", "run", "base", "contributing"])
    saved = saved[saved["fold"] == fold]
    agg, checks = None, []
    with tempfile.TemporaryDirectory(prefix="fig1b_dir_") as tmp:
        xp, yp, op = (os.path.join(tmp, n) for n in ("X.tsv", "y.tsv", "o.tsv"))
        X.to_csv(xp, sep="\t")
        y.rename("label").to_csv(yp, sep="\t", header=True)
        for run, seed in enumerate(SEEDS):
            r = subprocess.run(["Rscript", os.path.join(HERE, "all_directed.R"), xp, yp, str(seed), op],
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise RuntimeError(f"{dt} fold {fold} seed {seed}: {r.stderr[-400:]}")
            d = pd.read_csv(op, sep="\t")
            for c in ("base", "contributing"):
                d[c] = d[c].astype(str).str.replace(" ", "_")
            s = d[d["status"] == 2]
            score = {}
            for b, c_, t in zip(s["base"], s["contributing"], s["total"]):
                k = (b, c_) if b < c_ else (c_, b)
                if t > score.get(k, -1.0):
                    score[k] = t
            ranked = sorted(score.items(), key=lambda kv: -kv[1])
            top = {k for k, _ in ranked[:max(1, int(len(ranked) * 0.02))]}
            sv = saved[saved["run"] == run]
            ref = {tuple(sorted(p)) for p in zip(sv["base"], sv["contributing"])}
            checks.append({"data_type": dt, "fold": fold, "seed": seed, "n_directed_pairs": len(d),
                           "n_significant": int((d["status"] >= 1).sum()), "n_synergistic": len(s),
                           "top2pct_matches_rerun": top == ref})
            d["n"] = 1
            d["n_sig"] = (d["status"] >= 1).astype(int)
            d["n_syn"] = (d["status"] == 2).astype(int)
            part = d.groupby(["base", "contributing"])[["base_1D", "total", "n", "n_sig", "n_syn"]].sum()
            agg = part if agg is None else agg.add(part, fill_value=0)
    agg.reset_index().to_csv(os.path.join(OUT, f"{dt}_fold{fold:02d}.tsv.gz"), sep="\t", index=False,
                             compression="gzip", float_format="%.4f")
    return checks


def main():
    os.makedirs(OUT, exist_ok=True)
    only = sys.argv[1:]
    tasks = []
    for dt in ("taxa", "functions", "combined"):
        f = pd.read_csv(os.path.join(RESULTS, dt, f"mdfs_ig_info_{dt}.csv"),
                        usecols=["fold", "test_cohort"]).drop_duplicates().sort_values("fold")
        tasks += [(dt, int(a), b) for a, b in zip(f["fold"], f["test_cohort"]) if not only or f"{dt}:{a}" in only]
    checks = []
    with ProcessPoolExecutor(6) as ex:
        futs = {ex.submit(run_fold, *t): t for t in tasks}
        for i, fu in enumerate(as_completed(futs), 1):
            c = fu.result(); checks += c; t = futs[fu]
            print(f"[{i}/{len(tasks)}] {t[0]} fold {t[1]}: matches rerun in "
                  f"{sum(x['top2pct_matches_rerun'] for x in c)}/{len(c)} runs", flush=True)
    pd.DataFrame(checks).to_csv(os.path.join(HERE, "check_all_directed_test.tsv" if only else "check_all_directed.tsv"),
                                sep="\t", index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
