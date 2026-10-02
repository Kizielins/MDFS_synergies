# MDFS synergistic feature pairs in CRC and other microbiome cohorts — manuscript analysis code

This branch (`manuscript-analysis`) contains the exact code and result tables behind every figure,
table and number of the manuscript *"Synergy-oriented feature selection …"* (Zielińska, Rudnicki, Łabaj;
npj Biofilms and Microbiomes, revision). The stand-alone, single-cohort version of the method is the
tool on the [`main`](https://github.com/Kizielins/MDFS_synergies) branch.

Folders are named after the manuscript sections they reproduce (`M` = Materials and Methods,
`R` = Results, in manuscript order). Each notebook starts with a header that says which figures,
tables and quoted numbers it produces, and every manuscript item is written to `manuscript/`.

## Layout

| Folder | Manuscript section | Code | Produces |
|---|---|---|---|
| `run_mdfs.R` | Methods *MDFS feature selection and information gain* | one MDFS run: 1D IG + p-values, 2D significance, Bonferroni pair threshold, synergy filter | pair table, 1D table |
| `M1_data_preparation/` | Methods *Data preparation* | `prep_data.ipynb` | `prepared_data/X_*, y_*` (17 CRC cohorts × taxa / functions / combined; git-ignored) |
| `M2_loco_crc_pipeline/` | Methods *Leave-one-cohort-out (LOCO) validation* | `framework.ipynb` (functions), `main_code.ipynb` (orchestration), `run_loco.py` (headless) | `results/<dt>/`: `loco_results` (AUC), `loco_accuracy`, `base_features`, `synthetic_features` (RF selection frequencies), `mdfs_ig_info` (IGs of the top-2 % pairs of every fold × run), `mdfs_1d` (1D IG + p-values of every fold × run) |
| `R1_synergistic_pairs/` | Results *Synergy-oriented framework identifies disease-informative feature pairs* | `compare_mdfs_igs.ipynb`; `figure1b_all_pairs/` (all directed pairs, every fold × run); `seed_stability/` (20-seed stability grid) | **Figure 1b**, **Supp. Table 1**, **Supp. Table 5**, Wilcoxon test, `ig_gain_summary_<dt>.csv` |
| `R2_classification_performance/` | Results *Synergy-derived features retain benchmark-level classification performance* | `summarize_rf_results.ipynb`; `stable_pairs.ipynb`; `effect_sizes.ipynb` | AUC summary (`loco_performance_summary.tsv`), **Supp. Figures 1–2**, stable-pair table (`final_synthetic_features_table.csv`), **Figure 3**, **Supp. Table 2**, meta-analysis tables |
| `R3_pair_transformations/` | Results *Pair transformations reveal distinct forms of disease-associated dependence* | `transformation_preference.ipynb` | **Figure 4**, **Supp. Table 3**, LR / GM counts |
| `R4_broader_application/` | Results *Broader application identifies conserved synergistic markers in other conditions*; Methods *Application to additional disease cohorts* | `run_mdfs_1d_2d.py` (10 MDFS runs per cohort), `load_cohort_data.py` (synergy table), `analyse_synergy_gains.ipynb` | **Figure 5**, **Table 1**, **Supp. Table 4**, `synergy_by_cohort.tsv`, emergent synergies, top markers per cohort |
| `manuscript/` | — | — | one file per manuscript figure / table (written by the notebooks above, never by hand), `REVISION_CHANGES.md` (submitted vs revised numbers) |
| `archive/` | — | analyses not used in the manuscript (see `archive/README.md`) | — |
| `METHODS.md` | — | step-by-step description of the whole pipeline with all parameters | — |

## Method in one paragraph

For every LOCO fold (14 test-eligible of 17 CRC cohorts) MDFS is run **10 times** (seeds 42–51) on the
pooled training cohorts; in each run the pairs passing the Bonferroni pair-level χ² threshold and the
synergy filter (`total_IG > max(1D IGs)`) are ranked by total IG and the **top 2 %** kept; a pair is a
**consensus pair if it is in the top 2 % of ≥ 6 of the 10 runs**. Consensus pairs are encoded as
log-ratio (LR) and geometric-mean (GM) synthetic features and one Random Forest (`random_state = 42`)
is trained per feature configuration (raw / synthetic / raw + synthetic) and evaluated on the held-out
cohort. The validation arm applies the same MDFS filters to each non-CRC cohort (10 runs, seeds 12–21,
pair kept if it passes in ≥ 6 of 10 runs). Details and all parameters: `METHODS.md`.

## Reproducing everything

```bash
bash run_all.sh            # all steps in manuscript order; or e.g.  bash run_all.sh R1 R2 R3 R4
```

Steps: `M1` data preparation (needs the benchmark profiles in `../data/`), `M2` LOCO pipeline
(≈ 2–3 h, the three data types run in parallel), `R1`–`R4` the result notebooks (minutes each, except
`R2` effect sizes ≈ 1 h). Notebooks are executed in place with `jupyter nbconvert --execute`, so the
stored outputs are those of the last full run. The all-pairs input of Figure 1b
(`R1_synergistic_pairs/figure1b_all_pairs/run_all_directed.py`, several hours) is not part of
`run_all.sh`; its aggregated per-fold tables are deposited with the data (below).

## Data and large files (not in git)

- `M1_data_preparation/prepared_data/` and `R4_broader_application/data/`: input profiles
  (CRC benchmark of Piccinno et al. 2025 and curatedMetagenomicData) — see Data Availability in the paper.
- Per-run MDFS tables (`mdfs_ig_info_*`, `mdfs_1d_*`, `figure1b_all_pairs/directed/`): 15–95 MB each,
  regenerated by `M2`/`R4` or available from the data deposit. All summary tables and figures are in git.

## Software

R 4.5 with `MDFS` 1.5.5 and `data.table`; Python 3.9+ with `pandas`, `numpy`, `scikit-learn`, `scipy`,
`statsmodels`, `matplotlib`, `seaborn`, `jupyter` (`nbconvert`, `ipykernel`). Versions used for the
manuscript are listed in `METHODS.md` §7.
