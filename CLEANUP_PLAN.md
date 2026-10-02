> Status (2026-10-02): executed. Layout, path edits and notebook re-execution done; see README.md. Open decisions are listed in the chat summary.

# Repository clean-up plan (manuscript-analysis branch)

Drafted 2026-10-02. Nothing has been modified yet. Manuscript reference: `../../../NPJ_review/MS_final_reviewed.docx`.

## 0. Where we stand

| Item | Status |
|---|---|
| MDFS 10 runs, discovery arm | Done in code (`main_code.ipynb`: `N_RUNS_PER_FOLD = 10`, seeds 42–51, top 2 %, consensus ≥ 6/10) and in the `results_*_v2` outputs (`pair_ig_count` out of 140 = 14 folds × 10 runs). Code change is **uncommitted**. |
| MDFS 10 runs, validation arm | Done (`run_mdfs_1d_2d.py`: 10 runs, seeds 12–21, kept if ≥ 6/10). Outputs in `apply_framework/results_v2`. Uncommitted. |
| Random Forest | Trained **once** per fold and scenario (`random_state = 42`). The per-run RF of the submitted version was removed because its output was never used. If "10 iterations of RF" means 10 RF seeds per fold with the AUC averaged, that is a new change (see decision D2). |
| Downstream notebooks | Still hard-wired to the old folders (`results_taxa/`, `TOTAL_FOLD_RUNS = 28`, markdown saying "3 runs per fold"); the v2 outputs were produced by `rerun_v2/run_downstream.py`, which rewrites paths on the fly. Embedded notebook outputs are stale (main_code.ipynb carries 17 MB of 3-run output). |
| Manuscript | Results §1–2 numbers updated to the 10-run results; Methods (LOCO, RF, validation-arm paragraphs) and Results §3–4 still describe the submitted version (2 runs / 3 seeds, 62 pairs, 133,747 pairs, Gemmiger). Supplementary tables xlsx and Supplementary Figures docx are the submitted versions. |

## 1. Target layout (folders named after manuscript sections)

```
framework_repository/
├── README.md                     map: manuscript section → folder → figure/table; how to run end to end
├── METHODS.md                    (= methods.md, updated for the final layout)
├── run_mdfs.R                    Methods "MDFS feature selection and information gain"
├── requirements.txt, R_session_info.txt
├── run_all.sh                    runs every step below in order (one logged command per step)
│
├── M1_data_preparation/          Methods "Data preparation"
│   ├── prep_data.ipynb
│   └── prepared_data/            (git-ignored, unchanged)
│
├── M2_loco_crc_pipeline/         Methods "Leave-one-cohort-out (LOCO) validation"
│   ├── framework.ipynb           ML functions (unchanged logic)
│   ├── main_code.ipynb           orchestration, outputs cleared, paths fixed
│   ├── run_loco.py               headless runner (from rerun_v2/)
│   └── results/{taxa,functions,combined}/   loco_results, loco_accuracy, mdfs_ig_info,
│                                 base_features, synthetic_features, mdfs_1d_pvalues (new)
│
├── R1_synergistic_pairs/         Results "Synergy-oriented framework identifies disease-informative feature pairs"
│   ├── compare_mdfs_igs.ipynb    → Figure 1b, Supplementary Table 1 (ig_gain_summary per data type)
│   └── seed_stability/           → Supplementary Table 5 (scripts + consensus_grid.tsv copied from ../seed_stability_analysis)
│
├── R2_classification_performance/  Results "Synergy-derived features retain benchmark-level classification performance"
│   ├── summarize_rf_results.ipynb  → AUC/accuracy summary (Results text), Supp. Figure 1 (a,b), Supp. Figure 2
│   ├── stable_pairs.ipynb          (= consistent_mdfs_markers.ipynb minus the Fig 4 / Supp Table 3 cells) → stable-pair table (211)
│   └── effect_sizes.ipynb          (= consistent_mdfs_markers_downstream.ipynb) → Figure 3, Supplementary Table 2
│
├── R3_pair_transformations/      Results "Pair transformations reveal distinct forms of disease-associated dependence"
│   └── transformation_preference.ipynb  (Fig 4 + Supp Table 3 cells moved out of consistent_mdfs_markers) → Figure 4, Supplementary Table 3
│
├── R4_broader_application/       Results "Broader application identifies conserved synergistic markers in other conditions"
│   ├── data/                     (git-ignored, unchanged)
│   ├── run_mdfs_1d_2d.py, load_cohort_data.py
│   ├── results/                  (= apply_framework/results_v2)
│   ├── synergy_by_cohort.tsv     (= synergy_by_cohort_v2.tsv)
│   └── analyse_synergy_gains.ipynb → Figure 5, Table 1 (cohorts run / skipped), Supplementary Table 4, all counts quoted in Results §4
│
└── manuscript/                   one file per manuscript item, written by the notebooks above (never by hand)
    ├── Figure1a_schematic.pptx (= figure1.pptx), Figure1b.{png,pdf}, Figure3.*, Figure4.*, Figure5.*
    ├── SuppFigure1.*, SuppFigure2.*
    ├── SuppTable1..5.tsv  +  Supplementary_tables.xlsx (assembled from the TSVs)
    ├── RESULTS_NUMBERS.md        every number quoted in Results/Methods, with the producing notebook + cell
    └── REVISION_CHANGES.md       (= rerun_v2/REVISION_RESULTS.md, with the two known wrong rows corrected)
```

Naming: `M<n>` = Methods subsection order, `R<n>` = Results subsection order in the manuscript.

## 2. What gets deleted (git history keeps everything that was committed)

| Path | Why |
|---|---|
| `results_taxa/`, `results_functions/`, `results_combined/` | submitted 2-run results, superseded by `*_v2` |
| `apply_framework/results/`, `synergy_by_cohort.tsv`, `emergent_synergies_table.tsv`, `figure5_*`, `top_markers_per_cohort.txt`, `universal_disease_markers.txt`, `*.run.log`, `__pycache__/` | submitted 3-run validation-arm results and copies of its outputs |
| `downstream_analyses/` (outputs only; notebooks move) and `downstream_analyses_v2/` | superseded; regenerated into R2/R3 |
| `rerun_v2/loco_*.log` (103 MB), `work_*`, `test_original/`, `outputs/` | scaffolding of the September rerun; regenerated by the final run |
| `rerun_v2/figure1b_all_pairs/` (481 MB) | keep only the two scripts + the final figure if the "all pairs" Figure 1b is adopted (decision D4); otherwise delete |
| `run_taxa_analysis.py` | stub that does nothing |
| `RUN_TAXA_ANALYSIS.md` | describes the 3-run version; replaced by README |
| `compute_synergies.py`, `test_files/` | the stand-alone tool; lives on `main`, stale copy here |
| `investigate_synergies.r`, `pairs.csv`, `ig_synergy_*.csv` | hand-curated side check, not used by any manuscript figure/table, no seed (decision D6) |
| `analyse_synthetic_combinations.ipynb` (4.7 MB) | "incomplete pairings" forest plots; not used by any manuscript figure/table — the xlsx Supp Table 3 comes from consistent_mdfs_markers cell 6 (decision D6) |
| embedded outputs of `main_code.ipynb` (17 MB) | stale 3-run output |

## 3. Code changes

1. `run_mdfs.R` — (if D1 = yes) use one seed and the 2D discretisation range for the 1D IGs, the 2D significance test and `ComputeInterestingTuples`, as in the tool on `main`; keep the 5-output interface and keep writing 1D p-values (`MDFS(dimensions = 1, seed, range)`).
2. `framework.ipynb` — keep the 1D p-values of every run (`mdfs_1d_pvalues_<dt>.csv`) instead of deleting the temp file (needed for the P. micra 1D-significance statement, Reviewer 2 comment 4).
3. `main_code.ipynb` — write to `M2_loco_crc_pipeline/results/<dt>/`; run all three data types; clear outputs; (if D2 = yes) 10 RF seeds per fold, per-seed and mean AUC saved.
4. All downstream notebooks — read/write the final folders directly so `run_downstream.py`'s path-substitution is no longer needed; replace "3 runs"/"2 runs"/`TOTAL_FOLD_RUNS = 28` with the real values; each manuscript figure/table saved by an explicitly labelled cell into `manuscript/`.
5. `analyse_synergy_gains.ipynb` — count unordered pairs (max-total-IG orientation), as agreed in `NPJ_review/outstanding_mismatches.md`; Figure 5, the Results §4 numbers and Supp Table 4 all from the same table; emit the Table 1 run/skip status per cohort.
6. `stable_pairs` / `effect_sizes` — add the explicit Supp Table 2 cell (the old TSV was pasted via clipboard); report the GM/LR split under the stated definition.
7. New: `manuscript/RESULTS_NUMBERS.md` generator, `run_all.sh`, README, `.gitignore` additions (`*.log`, temp files).
8. Supplementary Table 5 from `seed_stability/consensus_grid.tsv`.

## 4. Execution order (after decisions)

1. Restructure (git mv), fix paths, delete superseded folders.
2. Run `prep_data` check (no changes expected) → LOCO ×3 data types (≈ 2–3 h in parallel) → validation arm (≈ 10 min).
3. Run downstream notebooks in place (`jupyter nbconvert --execute --inplace`) so the stored outputs are the final ones.
4. Build `manuscript/`, cross-check every number in `RESULTS_NUMBERS.md` against the docx, list the text edits still required.
5. Commit on `manuscript-analysis`; push after you have looked at the diff.

## 5. Issues found (details in the chat summary)

D1 discretisation mismatch in `run_mdfs.R` · D2 RF repeats · D3 Supplementary Table 4 empty under ordered grouping · D4 which Figure 1b · D5 git size of per-run IG tables · D6 what to archive · plus text items: Methods still says 2 runs / 3 seeds; Results §3 "41/62 (66 %)" and "FDR across the 62 pairs" stale; Results §4 counts stale; Figure 2 missing from numbering; Supp Table 5 cited but not yet produced; GM-share "127/211" in REVISION_RESULTS.md does not reproduce (161/211 under the submitted definition).
