#!/usr/bin/env python3
"""Supplementary Table 5: MDFS seed-stability grid (Reviewer 1, comment 6).

Reads consensus_grid.tsv (written by analyse_seed_stability.py from the 20-seed runs of
run_seed_stability.py) and writes ../../manuscript/SuppTable5_seed_stability.tsv:
for each data type, number of MDFS runs N and consensus threshold k, the mean consensus-set
size and the mean / minimum Jaccard agreement between consensus sets built from two disjoint
blocks of seeds (42-51 vs 52-61), averaged over the 14 LOCO folds.
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
g = pd.read_csv(os.path.join(HERE, "consensus_grid.tsv"), sep="\t")
summ = (g.groupby(["data_type", "N", "k"], sort=False)
          .agg(consensus_pairs_per_fold=("size", "mean"), jaccard_mean=("jaccard_AB", "mean"),
               jaccard_min=("jaccard_AB", "min"), n_folds=("fold", "nunique"))
          .reset_index())
summ["consensus_pairs_per_fold"] = summ["consensus_pairs_per_fold"].round(0).astype(int)
summ[["jaccard_mean", "jaccard_min"]] = summ[["jaccard_mean", "jaccard_min"]].round(3)
single = g.drop_duplicates(["data_type", "fold"]).groupby("data_type")["single_run_jaccard"].mean().round(3)
out = os.path.join(HERE, "..", "..", "manuscript", "SuppTable5_seed_stability.tsv")
summ.to_csv(out, sep="\t", index=False)
print(summ.to_string(index=False))
print("\nmean Jaccard between two single runs:\n", single.to_string())
print("\nSaved", os.path.abspath(out))
