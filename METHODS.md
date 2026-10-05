# Methods

This document records, step by step, how the CRC pooling analysis was conducted in this repository, with enough detail to reproduce or revise the pipeline. It describes the **revision** version of the code (10 MDFS runs per fold/cohort, top 2 %, consensus ≥ 6/10). Folder names follow the manuscript sections (`M` = Methods, `R` = Results); see `README.md`. What changed relative to the submitted version is summarised in §0 and in `manuscript/REVISION_CHANGES.md`.

The analysis has two independent arms:

1. **Discovery arm** — leave-one-cohort-out (LOCO) cross-validation on 17 colorectal cancer (CRC) cohorts to evaluate whether MDFS-derived synthetic feature pairs improve a Random Forest classifier (`M2_loco_crc_pipeline/framework.ipynb` + `main_code.ipynb`, run for each data type with `M2_loco_crc_pipeline/run_loco.py`).
2. **Validation arm** — applying the MDFS recipe per cohort to non-CRC microbiome cohorts to surface synergistic feature pairs across diseases (`R4_broader_application/`).

Both arms share `run_mdfs.R` for MDFS feature selection. Downstream notebooks (`R1_synergistic_pairs/compare_mdfs_igs.ipynb`, `R2_classification_performance/summarize_rf_results.ipynb`, `R2_classification_performance/stable_pairs.ipynb`, `effect_sizes.ipynb`, `R3_pair_transformations/transformation_preference.ipynb`, `R4_broader_application/analyse_synergy_gains.ipynb`) consume the outputs; `run_all.sh` runs everything in manuscript order (§5).

A stand-alone, single-cohort version of the same MDFS + RF procedure is on the `main` branch (`run_analysis.py`).

---

## 0. Changes from the submitted version

| | Submitted | Revision (this document) |
|---|---|---|
| MDFS runs per LOCO fold | 2; the seed was **not passed** to `run_mdfs.R`, so both runs were identical | **10 independent runs**, seeds 42–51 |
| Pairs kept per run | top 1 % by `total_IG` | **top 2 %** by `total_IG` |
| Consensus (LOCO) | pair in both runs | pair in **≥ 6 of 10** runs |
| Per-run diagnostic RF | trained each run, output (`rf_features_per_fold_<dt>.csv`) not used for any result | removed |
| Final RF | `random_state = 42` | unchanged |
| Validation arm | 3 runs (seeds 12–14), pair reported if significant in **any** run | **10 runs** (seeds 12–21), pair reported if significant in **≥ 6 of 10** runs |
| Feature / sample order | from Python sets (varied between sessions) | sorted (reproducible) |

The choice of 10 runs, top 2 % and ≥ 6/10 is based on the seed-stability analysis in `R1_synergistic_pairs/seed_stability/` (`consensus_grid.tsv`).

---

## 1. Input data and labels

`M1_data_preparation/prep_data.ipynb` produces, for each CRC cohort and each data type (`taxa`, `functions`, `combined`), the following files in `M1_data_preparation/prepared_data/`:

- `X_<cohort>_<data_type>.tsv` — samples × features matrix, tab-separated, first column = sample ID.
- `y_<cohort>_<data_type>.tsv` — sample ID + binary label, where **0 = CRC / disease**, **1 = control / healthy**.

The 17 cohort identifiers are listed as `ALL_COHORT_STEMS` in `main_code.ipynb`. They span:

- 11 public CRC studies: FengQ_2015, GuptaA_2019, LiuNN_2022, ThomasAM_2018b, VogtmannE_2016, WirbelJ_2018, YachidaS_2019, YangJ_2020, YangY_2021, YuJ_2015, ZellerG_2014.
- 6 in-house cohorts: ONCOBIOME_AtezoTRIBE, ONCOBIOME_COLOBIOME, ONCOBIOME_IIGM_CZ, ONCOBIOME_IIGM_IT, NSHII, IIGM_TU.

**Canonical cohort counts.** Different analyses apply different eligibility rules, so three denominators appear in this document and in the results files. They are not inconsistent — quote the one that matches the analysis being described:

| N | Meaning | Where it applies |
|---|---|---|
| **17** | all CRC cohorts loaded into `prepared_data/` | §1, and training sets throughout (each LOCO fold trains on the other 16) |
| **15** | cohorts with both label classes present | §5.4 cross-cohort meta-analysis (`n_cohorts` in `meta_results_all_pairs.csv`) |
| **14** | cohorts with ≥15 CRC **and** ≥15 controls | §3.1 LOCO test cohorts / folds (rows of `loco_results_<dt>.csv`) |

Separately, **140** is not a cohort count: it is 14 folds × 10 MDFS runs, the denominator of `pair_ig_count` in `ig_gain_summary_<dt>.csv` (§5.1).

Three data types per cohort:

- **`taxa`** — relative abundances of bacterial species (1,300–2,000 features per cohort).
- **`functions`** — pathway abundances (~500 features per cohort).
- **`combined`** — taxa and functions concatenated (~1,800–2,500 features per cohort).

For the **validation arm** (`R4_broader_application/data/`), cohorts are loaded from `merged_metadata.tsv`, `merged_abundances.txt`, `merged_pathways.txt`. Labels are derived from `study_condition` / `disease`: control/healthy → 1, any other condition → 0. CRC cohorts used by the discovery arm are excluded (see §4.1).

---

## 2. MDFS feature selection (`run_mdfs.R`)

For a given training matrix X and binary label vector y, `run_mdfs.R` runs the MDFS R package (v1.5.5) to identify (a) individually significant features and (b) pairs of features with significant **synergistic** information gain (IG) above and beyond the marginal 1D IG.

### 2.1 1D significance

```r
result_1d <- MDFS(X, y, dimensions = 1)
```

Default settings: `divisions = 1` (each variable is split in two bins), `discretizations = 1`, `pseudo.count = 0.25`, `p.adjust.method = "holm"`, `level = 0.05`. The set of variables flagged as 1D-relevant by MDFS is recorded along with each variable's IG, raw p-value, and Holm-adjusted p-value (written to `<output_1d_ig_file>` with columns `Feature, IG, p_value, adjusted_p_value`).

### 2.2 2D significance

```r
result_2d <- MDFS(X, y, dimensions = 2)
sig_2d_vars <- result_2d$relevant.variables
```

`MDFS(dimensions = 2)` flags variables for which there exists *at least one* partner producing a significant 2D IG at α = 0.05 after Holm correction across variables. This is a per-variable filter, not a per-pair filter.

### 2.3 Pair enumeration with chi²-based pair-level significance

The set of significant 2D variables is fed into `ComputeInterestingTuples`, which enumerates all pairs in which at least one endpoint is in `sig_2d_vars`:

```r
tuples <- ComputeInterestingTuples(
  X, y,
  dimensions   = 2,
  interesting.vars = sig_2d_vars,
  I.lower      = ig_1d,            # so reported tuples$IG = additional 2D IG over base's 1D
  ig.thr       = ig_thr_pair       # see derivation below
)
```

Setting `I.lower = ig_1d` means the IG returned for a tuple is the *additional* IG the contributing variable provides over the base variable's marginal 1D IG (i.e. the synergy contribution).

**Pair-level chi-squared cutoff (Bonferroni-corrected over the full pair search space):**

`MDFS::ComputePValue` (which we read in source) shows that under H₀, for default 2D MDFS parameters, the test statistic `2·log(2)·IG` follows χ² with

```
df = response.divisions × divisions × (divisions + 1)^(dimensions − 1)
   = 1 × 1 × 2^(2−1)
   = 2
```

We applied a Bonferroni correction over all theoretically possible unordered pairs in the matrix:

```r
n_total      <- ncol(X)
n_all_pairs  <- n_total × (n_total − 1) / 2
alpha        <- 0.05
ig_thr_pair  <- qchisq(1 − alpha / n_all_pairs, df = 2) / (2 · log(2))
```

Pairs with `IG_2D_added < ig_thr_pair` are dropped at this stage. (`run_mdfs.R` accepts an optional 7th argument `no_bonferroni` that skips this cutoff; neither arm uses it.)

### 2.4 Symmetric synergy filter

Pair-level significance can still be passed by pairs of the form (one strongly informative variable, one noise variable) because the joint contingency easily clears the threshold once the strong variable enters. To require **true synergy** — meaning the pair carries more information than either variable alone — we add:

```r
out_df <- out_df[out_df$total_IG > pmax(out_df$base_IG_1D, out_df$contributing_IG_1D), ]
```

The output written to `<output_csv>` has one row per (base, contributing) tuple with columns:

```
base, contributing, base_IG_1D, contributing_IG_1D, IG_2D_added, total_IG
```

where `total_IG = base_IG_1D + IG_2D_added`. Each unordered pair appears in **both** directions.

The unique 1D + 2D significant features list is also written to `<output_unique>`.

### 2.5 Seeds and repeated runs

MDFS discretisation has a random component, so a pair can pass the filters in some runs and not in others. `run_mdfs.R` takes an optional 6th positional argument, a random seed, used in `set.seed(seed)` at the top of the script (default 12 if omitted); MDFS draws its discretisation seed from R's random number generator, so a given seed gives identical results. Both arms run MDFS **10 times** with consecutive seeds and build a consensus:

- **LOCO arm** — seeds **42–51** (`seed = 42 + run_idx`), pair kept if in the per-run top 2 % in ≥ 6 of 10 runs (§3.3–3.6).
- **Validation arm** — seeds **12–21** (`seed = 12 + run_idx`), pair kept if significant in ≥ 6 of 10 runs (§4.2).

Training features are sorted and training samples concatenated in a fixed order before each run, so results are reproducible across sessions.

---

## 3. CRC LOCO discovery pipeline (`main_code.ipynb`)

`main_code.ipynb` orchestrates LOCO cross-validation over the 17 CRC cohorts for one data type at a time, calling the functions in `framework.ipynb`. For the revision, each data type was run with

```bash
python M2_loco_crc_pipeline/run_loco.py taxa        # likewise: functions, combined  (or: bash run_all.sh M2)
```

which executes the code cells of `framework.ipynb` and `main_code.ipynb` (except the notebook-loading and run-all cells) in a separate working directory per data type (`M2_loco_crc_pipeline/work_<dt>/`, so MDFS temporary files of parallel runs cannot collide) and writes the results to `M2_loco_crc_pipeline/results/<dt>/`. Logs: `M2_loco_crc_pipeline/loco_<dt>.log`.

### 3.1 Eligibility

Eligibility governs **which cohorts may be held out as a test set**. It does *not* restrict training: every fold trains on all 16 non-test cohorts, so all 17 cohorts — including the three listed below — enter the training data of at least 13 of the 14 folds.

A cohort is **test-eligible** iff its label vector contains at least 15 CRC samples and at least 15 controls. Of the 17 cohorts, **14 are test-eligible**:

| Cohort | CRC (0) | Control (1) | Excluded as test cohort because | In training folds? |
|---|---|---|---|---|
| ONCOBIOME_AtezoTRIBE | 163 | 0 | single-class — 0 controls, fails the ≥15-control rule | yes, all 14 |
| ONCOBIOME_COLOBIOME | 201 | 0 | single-class — 0 controls, fails the ≥15-control rule | yes, all 14 |
| NSHII | 14 | 448 | explicitly designated training-only via `TRAINING_ONLY_COHORT` | yes, all 14 |

NSHII is skipped by the explicit constant before the count rule is ever evaluated; it would also fail that rule independently (14 CRC < 15). The two ONCOBIOME cohorts are excluded purely by the count rule. Because none of the three is ever held out, each appears in the training set of **all 14 folds**.

This yields the **14 LOCO folds** listed in `loco_results_<data_type>.csv`. Note that this denominator is *not* the same as the **15** cohorts used in the cross-cohort meta-analysis of §5.4, which has a weaker eligibility requirement — see the note there before quoting either number.

### 3.1.1 Label–cohort coupling in the training set (known limitation)

Because the three cohorts above are retained in training, each training fold contains a subset of samples whose cohort of origin is almost perfectly predictive of the label:

- **AtezoTRIBE + COLOBIOME** contribute **364 samples, all CRC, no controls** — 12.6–15.2% of each training set, and **25.9–31.3% of all CRC samples in training**.
- **NSHII** contributes 448 controls and 14 CRC — **30.0–36.5% of all control samples in training**.

Together these three cohorts make up roughly **28–35% of each training set** while carrying near-degenerate label distributions. Any feature that separates these cohorts from the rest for technical or population reasons (sequencing batch, extraction protocol, geography) is therefore available to MDFS and the RF as an apparent CRC signal. This affects feature *discovery*; it does not leak into the reported AUCs, since every test cohort is a held-out cohort with both classes present and is never among the three. The relevant sensitivity analysis — refitting with the single-class cohorts dropped from training — is not currently in this repository.

### 3.2 Per-fold training set

For fold `i` testing on cohort `c`:

- training data = all 16 other cohorts, selected as `[s for s in ALL_COHORT_STEMS if s != test_cohort_stem]` — i.e. the test-eligibility rule of §3.1 is *not* applied here, so the single-class and training-only cohorts are included (see §3.1.1).
- features are intersected across the 16 training cohorts and sorted: 473–682 common features for taxa, 445–461 for functions, 918–1,143 for combined (range over the 14 folds).
- training matrices are concatenated, sample-wise.
- training set size = 2,393–2,884 samples.

The test cohort is reindexed to the training features (features missing in the test cohort are filled with 0).

### 3.3 Per-run MDFS + top-2 % filter

For each fold, MDFS is run **`N_RUNS_PER_FOLD = 10`** times on the merged training matrix, with seeds `42 + run_idx` (42–51). Each run:

1. Calls `run_mdfs.R` on the merged training matrix and labels (chi² + symmetric synergy filter, §2.3–2.4), via `run_mdfs_feature_selection` in `framework.ipynb`.
2. Deduplicates the returned tuples to unordered pairs and scores each pair by the larger `total_IG` of its two directions.
3. Keeps the **top 2 %** of pairs by that score (`TOP_PAIR_FRACTION = 0.02`, at least one pair).

Per run, this keeps (range over folds and runs):

| Data type | Pairs passing the MDFS filters | Kept (top 2 %) |
|---|---|---|
| taxa | 24,410–60,294 | 488–1,205 |
| functions | 46,384–59,473 | 927–1,189 |
| combined | 145,358–227,189 | 2,907–4,543 |

The IG values of the kept pairs (both directions) are recorded for every fold and run in `mdfs_ig_info_<dt>.csv` (§3.8).

### 3.4 Consensus across runs

After the 10 runs, a pair is a **consensus pair** if it was among the top 2 % in **at least 6 of the 10 runs** (`CONSENSUS_MIN_RUNS = 6`). Consensus pairs per fold (range, median over the 14 folds):

| Data type | Consensus pairs per fold |
|---|---|
| taxa | 533–1,104 (median 560) |
| functions | 326–692 (median 473) |
| combined | 3,087–4,494 (median 3,338) |

The consensus pairs are the only pairs used to build synthetic features for the final evaluation. No RF is trained during the 10 runs.

### 3.5 Synthetic feature generation

For each consensus pair `(f1, f2)` (names sorted alphabetically), two synthetic features are constructed (per sample, with `epsilon = 1e-9` to avoid division by zero):

- **`LR_<f1>__<f2>` = log((f1 + ε) / (f2 + ε))** — log ratio.
- **`GM_<f1>__<f2>` = sqrt((f1 + ε) · (f2 + ε))** — geometric mean.

So each fold adds twice as many synthetic features as it has consensus pairs (e.g. ~1,100 for taxa). The same pairs are used to compute the synthetic features of the test cohort.

### 3.6 Final evaluation — three scenarios

A Random Forest is trained from scratch on three feature spaces, each evaluated on the held-out test cohort:

| Scenario | Features |
|---|---|
| `All Features` | the base feature matrix only (no synthetic features) |
| `Synthetic Only` | the LR_/GM_ synthetic features derived from the consensus pairs |
| `All Features + Synthetic` | concatenation of the two |

Each scenario uses `train_and_evaluate()` (`framework.ipynb`), which does internal feature selection by:

1. Fitting an RF on a 75/25 stratified split of the training data (`random_state = 42`).
2. Ranking features by `feature_importances_`.
3. Validating subsets of sizes `k ∈ TOP_K_VALUES = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]` on the 25 % holdout, picking the `k*` with the highest validation AUC.
4. Ranking features again with an RF on the full training set and refitting the final RF on the top-`k*` features.

Test-set AUC and accuracy are recorded. A per-fold AUC summary is printed after each fold. The selected features of the `All Features + Synthetic` model are counted across folds (base and synthetic separately, §3.7).

### 3.7 Outputs (per data type)

Written to `M2_loco_crc_pipeline/results/<data_type>/`:

- `loco_results_<data_type>.csv` — AUC per test cohort × scenario, plus mean and weighted-by-N mean.
- `loco_accuracy_<data_type>.csv` — accuracy per test cohort × scenario.
- `mdfs_ig_info_<data_type>.csv` — IG values of the per-run top-2 % pairs for every (fold, run) in long format: `fold, test_cohort, run, base, contributing, base_IG_1D, contributing_IG_1D, IG_2D_added, total_IG` (`run` = 0–9).
- `base_features_<data_type>.csv` — how often each base feature was selected by the final RF in the *All Features + Synthetic* scenario across folds (`Feature, Selection Count, Selection %`).
- `synthetic_features_<data_type>.csv` — same for synthetic features.
- (`ig_gain_summary_<data_type>.csv` is written by `R1_synergistic_pairs/compare_mdfs_igs.ipynb` into `R1_synergistic_pairs/`, §5.1; `mdfs_1d_<data_type>.csv` holds the 1D IG, p-value and Holm-adjusted p-value of every feature for every fold and run.)

(`rf_features_per_fold_<data_type>.csv` of the submitted version is no longer produced.)

---

## 4. Validation arm — `R4_broader_application/`

The validation arm applies the MDFS recipe from §2 to non-CRC cohorts to ask: "do the same kinds of synergistic pairs emerge in other diseases / populations?" No classifier is trained in this arm.

### 4.1 Inputs

- `R4_broader_application/data/merged_metadata.tsv`, `merged_abundances.txt`, `merged_pathways.txt`.
- Cohorts are the unique `study_name` / `source_cohort` values in the metadata. The 11 CRC studies of the discovery arm plus HanniganGD_2017 and BaxterNT_2016 are excluded (`EXCLUDED_CRC_COHORTS` in `load_cohort_data.py`), leaving **51 cohorts**.
- For each cohort, taxa and pathways are concatenated over the samples present in both (sorted). A cohort is skipped if it has no data, fewer than 10 samples, or only one label class: **38 are skipped, 13 are analysed**.

### 4.2 Per-cohort MDFS (`run_mdfs_1d_2d.py`)

For each analysed cohort:

1. Run `run_mdfs.R` **10 times** (`N_RUNS_PER_COHORT = 10`) with seeds 12–21, with the Bonferroni pair-level cutoff and synergy filter (§2.3–2.4). No top-fraction cutoff is applied.
2. Keep (base, contributing) directions that pass the filters in **at least 6 of 10 runs** (`MIN_RUNS_FOR_PAIR = 6`, defined in `load_cohort_data.py`) and average their IGs over the runs in which they passed.
3. Write to `R4_broader_application/results/`:
   - `mdfs_ig_info_<cohort>.tsv` — all rows from all runs.
   - `ig_gain_summary_<cohort>.tsv` — per retained direction: `base_IG_1D_mean`, `contributing_IG_1D_mean`, `IG_2D_added_mean`, `IG_2D_added_std`, `total_IG_mean`, `n_runs`, `ig_gain` (= `IG_2D_added_mean`), `ig_gain_pct = IG_2D_added_mean / base_IG_1D_mean × 100`.
   - `mdfs_2d_<cohort>.tsv`, `mdfs_1d_ig_<cohort>.tsv` — raw output of the last run.

Of the 13 analysed cohorts, 7 have at least one retained pair (QinN_2014, HallAB_2017, JieZ_2017, MetaCardis_2020_a, RubelMA_2020, NielsenHB_2014, ZhuF_2020). Log: `R4_broader_application/validation_mdfs.log`.

```bash
python R4_broader_application/run_mdfs_1d_2d.py --data-dir R4_broader_application/data --all-cohorts --out-dir R4_broader_application/results
```

### 4.3 Cross-cohort synergy table (`load_cohort_data.py`)

`load_cohort_data.py` aggregates all `mdfs_ig_info_*.tsv` files into a single table grouped by `(cohort, base, contributing)`, applying the same ≥ 6/10 rule, with mean IG values and `ig_gain` / `ig_gain_pct` columns: `R4_broader_application/synergy_by_cohort.tsv` (80,110 cohort × direction rows).

```bash
python R4_broader_application/load_cohort_data.py --data-dir R4_broader_application/data --mdfs-results-dir R4_broader_application/results --out-table R4_broader_application/synergy_by_cohort.tsv
```

### 4.4 Cohort-level reporting (`analyse_synergy_gains.ipynb`)

Reads the synergy table and reduces it to **one row per cohort and unordered pair** (MDFS reports each pair in both orientations; the orientation with the larger `total_IG_mean` is kept and `pair_key` identifies the pair — all counts in Results §4 are unique pairs), plots base 1D IG vs total IG colored by `ig_gain_pct`, filters to `ig_gain_pct ≥ 50` for downstream, and writes:

- `manuscript/Figure5.{pdf,png}` — Figure 5.
- `manuscript/Table1_validation_cohorts.tsv` — Table 1: every cohort of the metadata with sample counts, whether it was analysed and how many unordered pairs it yielded.
- `top_markers_per_cohort.txt` — top 20 pairs per cohort by total_IG.
- `manuscript/SuppTable7_conserved_pairs.tsv` (= `universal_disease_markers.txt`) — unordered pairs found in ≥ 2 cohorts with the same disease label (Supp. Table 7).
- `emergent_synergies_table.tsv` — pairs in which one feature has zero 1D IG (quoted in Results §4).

---

## 5. Downstream notebooks

Each notebook reads the final result folders directly and is executed in place (`jupyter nbconvert --execute --inplace`, see `run_all.sh`), so the outputs stored in the notebooks are those of the last full run. Manuscript figures and tables are written to `manuscript/`.

### 5.1 `R1_synergistic_pairs/compare_mdfs_igs.ipynb`

Cross-data-type comparison of the LOCO IG outputs. Reads `M2_loco_crc_pipeline/results/<dt>/mdfs_ig_info_<dt>.csv` for taxa, functions, combined; aggregates per (base, contributing) with mean `base_IG_1D`, `total_IG` and their difference (= `IG_2D_added`); tests 2D vs 1D IG with a Wilcoxon signed-rank test; writes `R1_synergistic_pairs/ig_gain_summary_<dt>.csv` and `manuscript/SuppTable2_top30_synergistic_pairs.tsv`, and plots `base_IG_1D` vs `total_IG` of the consensus pairs (`figure1b_consensus_pairs_only.png`). The manuscript Figure 1b shows *all* directed pairs of every fold and run (significant or not, synergistic or not) and is produced by `R1_synergistic_pairs/figure1b_all_pairs/` (`run_all_directed.py` → `directed/<dt>_fold<NN>.tsv.gz`, `plot_figure1b_all_directed.py` → `manuscript/Figure1b.{png,pdf}`).

The `pair_ig_count` column of `ig_gain_summary_<dt>.csv` reports in how many fold–run combinations the (base, contributing) direction was among the per-run top 2 %, out of **140 = 14 LOCO folds × 10 MDFS runs per fold** (§3.3). A value of `140/140` therefore means the pair was in the top 2 % in every fold under every seed — the maximum possible consistency, not a cohort count.

### 5.2 `investigate_synergies.r`

Stand-alone R script that, for a hand-curated `pairs.csv` (a small list of pairs the analyst wants to scrutinize), computes 1D and 2D IGs of each listed pair on the merged training data of a chosen data type. Outputs `ig_synergy_<data_type>.csv` with columns `IG_1D_Feature1, IG_1D_Feature2, IG_2D_Pair, IG_Gain_vs_Max1D, IG_Gain_vs_Avg1D` plus 1D IG of the LR/MULT/GM transforms of each pair. Run for taxa, functions, and combined separately. It does not set a random seed and was not rerun for the revision.

`pairs.csv` is hand-curated — not automatically derived from the MDFS output.

### 5.3 `consistent_mdfs_markers.ipynb`

Reads `M2_loco_crc_pipeline/results/<dt>/{base_features,synthetic_features}_<dt>.csv` for all three data types. A synthetic pair is **stable** if

1. it was selected in ≥ 80 % of the 14 LOCO folds (`SYNTHETIC_THRESHOLD = 80`, i.e. ≥ 12 folds). The LR_ and GM_ features of a pair count as one pair, using the more frequently selected of the two (the first row after the input table's sort by selection count), and
2. its selection % is higher than that of each of its two component features as base features (which also makes the "selection % gain", pair % − mean of the two component %, positive).

Writes `R2_classification_performance/final_synthetic_features_table.csv`: **211 stable pairs** (combined 160, taxa 34, functions 17). It also writes Supplementary Table 6 (`supplementary_table_3_transformation_preference.tsv`: which transformation, LR or GM, reached ≥ 50 % selection; 746 pairs) and Figure 4 (`figure4_lr_vs_gm.{pdf,png}`).

### 5.4 `consistent_mdfs_markers_downstream.ipynb`

Takes the stable-pair table from §5.3, and for each pair computes:

- per-cohort Mann–Whitney U p-values (with BH-FDR correction within each pair across cohorts).
- per-cohort effect sizes (rank-biserial correlation, equivalent to AUC).
- across-cohort combined p-values (Stouffer's method, weighted by √N) with BH-FDR at the pair level.

**Cohort denominator — 15, not 14.** The notebook iterates over all 17 `ALL_COHORT_STEMS`, but a cohort contributes a test only if *both* label classes are present. AtezoTRIBE and COLOBIOME are single-class and drop out; NSHII does contribute here, because a Mann–Whitney test needs only ≥1 sample per class, not the ≥15-per-class LOCO threshold of §3.1. Hence:

```
17 cohorts
 − 2 single-class (AtezoTRIBE, COLOBIOME)   = 15  ← meta-analysis denominator (§5.4)
 − 1 under-powered (NSHII, 14 CRC < 15)     = 14  ← LOCO test-cohort denominator (§3.1)
```

This is confirmed by `R2_classification_performance/meta_results_all_pairs.csv`, where `n_cohorts = 15` for all 211 pairs. The `direction` column is reported as `+x/-y` with `x + y = 15`.

**Directional consistency.** Across the 211 pairs the observed direction splits are `+15/-0` (107 pairs), `+14/-1` (79), `+0/-15` (11), `+1/-14` (4), `+13/-2` (7), `+2/-13` (2) and `+12/-3` (1). So **201 of 211** pairs agree in direction in at least 14 of 15 cohorts (118 in all 15); 9 agree in 13 of 15 and 1 in 12 of 15. A blanket claim of "at least 14 of 15" applies to those 201 pairs only — state it as "at least 12 of 15 cohorts (14 or 15 for 201 of the 211 pairs)" if the claim is made over the full set.

Writes `meta_results_all_pairs.csv` and `meta_significant_pairs_fdr_lt_0.05.csv` (both contain the same 211 pairs — every pair tested passes FDR < 0.05), and Figure 3 (`figure3_hub_representative_pairs.png`). The overview figure with one panel per stable pair (`figure_synthetic_barplots.png`, cell 6) is not produced for the revision: with 211 pairs the image exceeds matplotlib's size limit.

---

## 6. Random Forest model details

Used identically for all RF training calls (3 scenarios × 14 folds × 3 data types):

```python
{
    'n_estimators'    : 1000,
    'max_features'    : 'sqrt',
    'min_samples_leaf': 1,
    'class_weight'    : 'balanced',
    'random_state'    : 42,
    'n_jobs'          : -1,
}
```

Internal feature-selection hyperparameter `k`:

- search grid: `[50, 100, 150, 200, 250, 300, 350, 400, 450, 500]`.
- selection criterion: AUC on a 25% stratified holdout (`random_state=42`).
- final RF refit on the top-`k*` features and evaluated on the held-out test cohort.

Performance metrics:

- AUC = `sklearn.metrics.roc_auc_score(y_true, y_pred_proba)` where `y_pred_proba = predict_proba(X_test)[:, 1]`.
- Accuracy is reported alongside AUC.
- Cohort-level results are aggregated as both unweighted mean and weighted (by test-cohort N) mean.

---

## 7. Software

- R 4.5.1, `MDFS` 1.5.5, `data.table` 1.18.4.
- Python 3.9.6, `pandas` 2.3.1, `numpy` 2.0.2, `scikit-learn` 1.6.1, `scipy` 1.13.1, `statsmodels` 0.14.5, `matplotlib` 3.9.4, `seaborn` 0.13.2, `IPython` 8.18.1 (used by `run_all.sh`).

---

## 8. Filter recipe — summary

Putting §2 and §3 together, a pair survives all the way through to the final RF evaluation only if:

1. At least one endpoint passes MDFS's per-variable 2D significance test (Holm-corrected α=0.05).
2. The pair's `IG_2D_added` exceeds the Bonferroni-corrected chi² threshold at α=0.05 over `n × (n−1) / 2` total possible pairs (with df=2, expressed in bits via `qchisq / (2·log(2))`).
3. The pair's `total_IG > max(base_1D_IG, contributing_1D_IG)` — true synergy.
4. The pair ranks in the top 2 % by `total_IG` within an MDFS run.
5. Condition 4 holds in at least 6 of the 10 MDFS runs of the fold (consensus).

Conditions (1)–(2) are statistical significance; (3) ensures the pair carries genuinely synergistic information; (4)–(5) bound the volume so the RF stage is tractable and keep only pairs that are reproducible across MDFS runs.

In the validation arm conditions (1)–(3) must hold in at least 6 of 10 runs; there is no top-2 % cutoff.

---

## 9. Reviewer checks (`R5_reviewer_checks/`)

Analyses added in the revision in response to the reviewers; every output is a TSV in the folder and the numbers are listed in `manuscript/RESULTS_NUMBERS.md` (section "R5 reviewer checks").

| Script / notebook | Reviewer comment | What it does | Output |
|---|---|---|---|
| `run_negative_controls.py` | R1.7 | LOCO fold-0 training set of each data type: `run_mdfs.R` with real labels, with labels permuted (5×), and with every feature shuffled independently within each class (5×; keeps univariate signal, destroys pair dependence) | `negative_controls.tsv` (Supp. Table 8) |
| `run_discretisation_sensitivity.py` | R1.10 | fold-0 taxa: default vs seed 43, wider (1D) cut range, half range, three bins (divisions = 2); overlap of pair sets with the default | `discretisation_sensitivity.tsv` (Supp. Table 9) |
| `run_runtime.py` | R1.9 | wall-clock time and peak memory of one MDFS run per data type (fold-0 training set) | `runtime.tsv` |
| `additivity_and_collinearity.ipynb` | R2.2, R2.5, R1.5 | interaction information `total_IG − base_IG_1D − contributing_IG_1D` of every per-run top-2 % pair and of the consensus pairs; cohort-adjusted logistic regression with likelihood-ratio tests (cohort-only vs additive vs interaction) for the stable pairs; Spearman correlation of each stable synthetic feature with its constituents | `interaction_information.tsv`, `stable_pairs_lrt.tsv` (Supp. Table 5), `collinearity.tsv` |
| `sensitivity_training_cohorts.ipynb` | R1.5 | compares the main LOCO run with `run_loco.py --exclude-training-only` (AtezoTRIBE, COLOBIOME and NSHII removed from every training set; results in `M2_loco_crc_pipeline/results_sensitivity/`): AUCs and stable pairs | `sensitivity_auc.tsv` (Supp. Table 3), `sensitivity_shared_stable_pairs.tsv` |
| `reviewer_checks_summary.ipynb` | — | tabulates the three script outputs | `r5_numbers_scripts.tsv` |

The 1D significance of the hub species (Results §2) is taken from the Holm-adjusted 1D MDFS p-values saved for every fold and run (`M2_loco_crc_pipeline/results/<dt>/mdfs_1d_<dt>.csv`, summarised by `R2_classification_performance/hub_1d_significance.py`).
