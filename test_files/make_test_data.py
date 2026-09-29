#!/usr/bin/env python3
"""
Generate the synthetic test dataset in test_files/ (X_test.tsv, y_test.tsv).

300 samples (150 disease = 0, 150 control = 1) x 30 bacterial species, as
relative abundances (each sample sums to 1).

Planted signal:
  - Synergistic pair: Parvimonas_micra and Gemella_morbillorum share a strongly
    varying sample-specific load, which blurs each one's individual signal,
    while their ratio P. micra / G. morbillorum is high in disease and low in
    controls. Each is somewhat informative alone, the pair much more so. MDFS
    should report this pair and its log-ratio (LR_) feature should improve the
    Random Forest.
  - Weak individual markers: Fusobacterium_nucleatum and
    Peptostreptococcus_anaerobius are slightly higher in disease,
    Faecalibacterium_prausnitzii and Roseburia_intestinalis slightly lower.
  - All other species are noise.

Run from the repository root:
    python test_files/make_test_data.py
"""

import os

import numpy as np
import pandas as pd

SEED = 2026
N_PER_CLASS = 150

SPECIES = [
    "Fusobacterium_nucleatum", "Peptostreptococcus_anaerobius", "Parvimonas_micra",
    "Gemella_morbillorum", "Clostridium_symbiosum", "Bacteroides_fragilis",
    "Prevotella_copri", "Faecalibacterium_prausnitzii", "Bifidobacterium_longum",
    "Ruminococcus_gnavus", "Akkermansia_muciniphila", "Roseburia_intestinalis",
    "Eubacterium_rectale", "Bacteroides_vulgatus", "Bacteroides_uniformis",
    "Alistipes_putredinis", "Escherichia_coli", "Streptococcus_salivarius",
    "Veillonella_parvula", "Dorea_longicatena", "Blautia_obeum",
    "Coprococcus_comes", "Anaerostipes_hadrus", "Collinsella_aerofaciens",
    "Parabacteroides_distasonis", "Bilophila_wadsworthia", "Odoribacter_splanchnicus",
    "Lachnospira_pectinoschiza", "Methanobrevibacter_smithii", "Dialister_invisus",
]

# Shift of the log-abundance in disease samples (weak individual markers)
MARKER_SHIFTS = {
    "Fusobacterium_nucleatum": 0.35,
    "Peptostreptococcus_anaerobius": 0.35,
    "Faecalibacterium_prausnitzii": -0.35,
    "Roseburia_intestinalis": -0.35,
}
PAIR = ("Parvimonas_micra", "Gemella_morbillorum")
PAIR_LOAD_SD = 3.0     # shared load: blurs each pair member's individual signal
PAIR_RATIO_SHIFT = 0.9  # +/- shift of each member's log-abundance, opposite directions


def main():
    rng = np.random.default_rng(SEED)
    n = 2 * N_PER_CLASS
    y = np.array([0] * N_PER_CLASS + [1] * N_PER_CLASS)
    disease = (y == 0)

    log_ab = pd.DataFrame(rng.normal(0.0, 1.0, size=(n, len(SPECIES))), columns=SPECIES)
    log_ab += rng.normal(0.0, 1.0, size=len(SPECIES))  # species-specific mean abundance

    for sp, shift in MARKER_SHIFTS.items():
        log_ab.loc[disease, sp] += shift

    a, b = PAIR
    load = rng.normal(0.0, PAIR_LOAD_SD, size=n)
    direction = np.where(disease, 1.0, -1.0)
    log_ab[a] += load + PAIR_RATIO_SHIFT * direction
    log_ab[b] += load - PAIR_RATIO_SHIFT * direction

    ab = np.exp(log_ab)
    rel = ab.div(ab.sum(axis=1), axis=0)

    sample_ids = [f"S{i + 1:03d}" for i in range(n)]
    rel.index = pd.Index(sample_ids, name="sample_id")
    labels = pd.Series(y, index=rel.index, name="label")

    out_dir = os.path.dirname(os.path.abspath(__file__))
    rel.to_csv(os.path.join(out_dir, "X_test.tsv"), sep="\t", float_format="%.8f")
    labels.to_csv(os.path.join(out_dir, "y_test.tsv"), sep="\t")
    print(f"Wrote {n} samples x {len(SPECIES)} features to {out_dir}")


if __name__ == "__main__":
    main()
