#!/usr/bin/env bash
# Reproduces every analysis of the manuscript from the prepared data, in manuscript order.
# Step 2 (LOCO) is the expensive part (about 2-3 h with the three data types in parallel on 6+ cores);
# every other step takes seconds to minutes, except R2 effect sizes (about 1 h).
# Usage: bash run_all.sh [step ...]      (no argument = all steps; e.g. "bash run_all.sh R1 R2 R3 R4")
set -euo pipefail
cd "$(dirname "$0")"
NB() { python3 -m jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=python3 "$1"; }
steps=("$@"); [ ${#steps[@]} -eq 0 ] && steps=(M1 M2 R1 R2 R3 R4)
has() { for s in "${steps[@]}"; do [ "$s" = "$1" ] && return 0; done; return 1; }

if has M1; then   # Methods: data preparation  (needs ../data/, the CRC benchmark profiles)
  NB M1_data_preparation/prep_data.ipynb
fi
if has M2; then   # Methods: LOCO validation  (10 MDFS runs per fold, one RF per configuration)
  for dt in taxa functions combined; do
    python3 M2_loco_crc_pipeline/run_loco.py $dt > M2_loco_crc_pipeline/loco_$dt.log 2>&1 &
  done; wait
fi
if has R1; then   # Results 1: synergistic pairs (Figure 1b, Supp Tables 1 and 5)
  NB R1_synergistic_pairs/compare_mdfs_igs.ipynb
  python3 R1_synergistic_pairs/seed_stability/make_supp_table5.py
  # Figure 1b with all pairs needs the per-run MDFS tables of every directed pair:
  #   python3 R1_synergistic_pairs/figure1b_all_pairs/run_all_directed.py   (hours; writes directed/*.tsv.gz)
  python3 R1_synergistic_pairs/figure1b_all_pairs/plot_figure1b_all_directed.py
fi
if has R2; then   # Results 2: RF performance, stable pairs, effect sizes (Supp Figs 1-2, Figure 3, Supp Table 2)
  NB R2_classification_performance/summarize_rf_results.ipynb
  NB R2_classification_performance/stable_pairs.ipynb
  NB R2_classification_performance/effect_sizes.ipynb
fi
if has R3; then   # Results 3: LR vs GM (Figure 4, Supp Table 3)
  NB R3_pair_transformations/transformation_preference.ipynb
fi
if has R4; then   # Results 4: other diseases (Figure 5, Table 1, Supp Table 4)
  ( cd R4_broader_application
    python3 run_mdfs_1d_2d.py --data-dir data --all-cohorts --out-dir results > validation_mdfs.log 2>&1
    python3 load_cohort_data.py --data-dir data --mdfs-results-dir results --out-table synergy_by_cohort.tsv )
  NB R4_broader_application/analyse_synergy_gains.ipynb
fi
echo "run_all.sh finished: ${steps[*]}"
