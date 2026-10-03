#!/usr/bin/env bash
# Overnight rerun with the fixed run_mdfs.R. Phases write markers to overnight_status.log.
cd "$(dirname "$0")"
REPO="/Users/kizie/Library/Mobile Documents/com~apple~CloudDocs/OneDrive-UJ/Sonata2/mdfs_combine_taxonomy_and_functions/CRC_segata_pooling/framework_repository"
cd "$REPO"
ST="$REPO/overnight_status.log"
mark() { echo "$(date '+%F %T')  $*" | tee -a "$ST"; }
NB() { python3 -m jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=python3 "$1"; }
export PYTHONUNBUFFERED=1
rm -f "$ST"; mark "START overnight run (branch $(git branch --show-current), $(git rev-parse --short HEAD))"

# ---------------- Phase A: LOCO x3 + validation arm ----------------
mark "A start: LOCO taxa/functions/combined + validation arm"
rm -rf M2_loco_crc_pipeline/results M2_loco_crc_pipeline/work_*
for dt in taxa functions combined; do
  ( python3 M2_loco_crc_pipeline/run_loco.py $dt > M2_loco_crc_pipeline/loco_$dt.log 2>&1; mark "A LOCO $dt exit=$?" ) &
done
( cd R4_broader_application && rm -rf results && python3 run_mdfs_1d_2d.py --data-dir data --all-cohorts --out-dir results > validation_mdfs.log 2>&1; mark "A validation MDFS exit=$?"
  python3 load_cohort_data.py --data-dir data --mdfs-results-dir results --out-table synergy_by_cohort.tsv > validation_table.log 2>&1; mark "A validation table exit=$?" ) &
wait
mark "A done"

# ---------------- Phase B: Figure 1b all pairs + sensitivity LOCO ----------------
mark "B start: figure1b all directed pairs (6 procs) + sensitivity LOCO (3 procs)"
rm -rf R1_synergistic_pairs/figure1b_all_pairs/directed M2_loco_crc_pipeline/results_sensitivity
( cd R1_synergistic_pairs/figure1b_all_pairs && python3 run_all_directed.py > run_all_directed.log 2>&1; mark "B figure1b all-directed exit=$?" ) &
for dt in taxa functions combined; do
  ( python3 M2_loco_crc_pipeline/run_loco.py $dt --exclude-training-only > M2_loco_crc_pipeline/loco_${dt}_sensitivity.log 2>&1; mark "B sensitivity LOCO $dt exit=$?" ) &
done
wait
mark "B done"

# ---------------- Phase C: seed stability (20 seeds x 42 folds) + reviewer checks ----------------
mark "C start: seed stability + reviewer checks"
( cd R1_synergistic_pairs/seed_stability && rm -rf pairs && python3 run_seed_stability.py > run.log 2>&1; mark "C seed stability exit=$?"
  python3 analyse_seed_stability.py > analyse.log 2>&1; mark "C seed analysis exit=$?"; python3 make_supp_table5.py > supp5.log 2>&1 ) &
if [ -d R5_reviewer_checks ]; then
  ( cd R5_reviewer_checks && for s in run_negative_controls.py run_discretisation_sensitivity.py run_runtime.py; do
      [ -f "$s" ] && { python3 "$s" > "${s%.py}.log" 2>&1; mark "C $s exit=$?"; }; done ) &
fi
wait
mark "C done"

# ---------------- Phase D: result notebooks, tables, numbers ----------------
mark "D start: notebooks"
NB R1_synergistic_pairs/compare_mdfs_igs.ipynb > R1_synergistic_pairs/nb.log 2>&1; mark "D R1 compare exit=$?"
python3 R1_synergistic_pairs/figure1b_all_pairs/plot_figure1b_all_directed.py > R1_synergistic_pairs/figure1b_all_pairs/plot.log 2>&1; mark "D fig1b plot exit=$?"
NB R2_classification_performance/summarize_rf_results.ipynb > R2_classification_performance/nb_summarize.log 2>&1; mark "D R2 summarize exit=$?"
NB R2_classification_performance/stable_pairs.ipynb > R2_classification_performance/nb_stable.log 2>&1; mark "D R2 stable exit=$?"
( NB R2_classification_performance/effect_sizes.ipynb > R2_classification_performance/nb_effect.log 2>&1; mark "D R2 effect exit=$?" ) &
( NB R3_pair_transformations/transformation_preference.ipynb > R3_pair_transformations/nb.log 2>&1; mark "D R3 exit=$?" ) &
( NB R4_broader_application/analyse_synergy_gains.ipynb > R4_broader_application/nb.log 2>&1; mark "D R4 exit=$?" ) &
wait
python3 R2_classification_performance/hub_1d_significance.py > R2_classification_performance/hub_1d.log 2>&1; mark "D hub 1D exit=$?"
if [ -d R5_reviewer_checks ]; then for nb in R5_reviewer_checks/*.ipynb; do NB "$nb" > "${nb%.ipynb}.log" 2>&1; mark "D $(basename $nb) exit=$?"; done; fi
python3 manuscript/collect_numbers.py > manuscript/collect_numbers.log 2>&1; mark "D numbers exit=$?"
mark "ALL DONE"
