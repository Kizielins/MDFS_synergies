"""Shared helpers for the reviewer-check scripts: LOCO training sets and one MDFS run."""
import os
import subprocess
import tempfile
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PREP = os.path.join(REPO, "M1_data_preparation", "prepared_data")
RUN_MDFS = os.path.join(REPO, "run_mdfs.R")
ALL_COHORT_STEMS = [
    'Public_study__FengQ_2015', 'Public_study__GuptaA_2019', 'Public_study__LiuNN_2022',
    'Public_study__ThomasAM_2018b', 'Public_study__VogtmannE_2016', 'Public_study__WirbelJ_2018',
    'Public_study__YachidaS_2019', 'Public_study__YangJ_2020', 'Public_study__YangY_2021',
    'Public_study__YuJ_2015', 'Public_study__ZellerG_2014', 'This_study_cohort1__ONCOBIOME_AtezoTRIBE',
    'This_study_cohort2__ONCOBIOME_COLOBIOME', 'This_study_cohort3__ONCOBIOME_IIGM_CZ',
    'This_study_cohort4__ONCOBIOME_IIGM_IT', 'This_study_cohort5__NSHII', 'This_study_cohort6__IIGM_TU'
]
FOLD0_TEST = 'Public_study__FengQ_2015'   # LOCO fold 0


def training_set(data_type, test_cohort=FOLD0_TEST):
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


def run_mdfs(X, y, seed, rscript=RUN_MDFS, extra_args=()):
    """One run of run_mdfs.R. Returns (pairs DataFrame, stdout text, wall seconds, max RSS bytes)."""
    import time
    with tempfile.TemporaryDirectory(prefix="r5_") as tmp:
        xp, yp = os.path.join(tmp, "X.tsv"), os.path.join(tmp, "y.tsv")
        X.to_csv(xp, sep="\t")
        y.rename("label").to_csv(yp, sep="\t", header=True)
        out = [os.path.join(tmp, n) for n in ("pairs.tsv", "uniq.txt", "ig1d.tsv")]
        t0 = time.time()
        r = subprocess.run(["/usr/bin/time", "-l", "Rscript", rscript, xp, yp, *out, str(seed), *extra_args],
                           capture_output=True, text=True)
        wall = time.time() - t0
        if r.returncode != 0:
            raise RuntimeError(r.stderr[-800:])
        rss = None
        for line in r.stderr.splitlines():
            if "maximum resident set size" in line:
                rss = int(line.split()[0])
        pairs = pd.read_csv(out[0], sep="\t") if os.path.getsize(out[0]) > 0 else pd.DataFrame(columns=["base", "contributing", "total_IG"])
        ig1d = pd.read_csv(out[2], sep="\t")
        return pairs, ig1d, r.stdout, wall, rss


def unordered(pairs):
    return {tuple(sorted((a, b))) for a, b in zip(pairs["base"], pairs["contributing"])} if len(pairs) else set()


def top_fraction(pairs, frac=0.02):
    if len(pairs) == 0:
        return set()
    d = pairs.copy()
    d["key"] = [tuple(sorted((a, b))) for a, b in zip(d["base"], d["contributing"])]
    sc = d.groupby("key")["total_IG"].max().sort_values(ascending=False)
    return set(sc.index[:max(1, int(len(sc) * frac))])
