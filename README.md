# Metagenomic Feature Synergy Analysis

Finds synergistic pairs of metagenomic features (taxa, pathways, or both) with
**Multidimensional Feature Selection (MDFS)**, turns the most robust pairs into
synthetic features, and tests whether they improve a Random Forest classifier.

The tool analyses a **single cohort**: one feature matrix and one set of
binary labels (disease vs control).

```bash
python run_analysis.py                                   # bundled test data
python run_analysis.py --X my_features.tsv --y my_labels.tsv --out-dir my_results
```

---

## What the tool does

```
                 ┌─────────────────────────────────────────────┐
 X, y ──────────►│ 1. MDFS on the full dataset (10 runs)       │──► mdfs_synergies.tsv
                 │    consensus pairs → LR_/GM_ features       │──► consensus_pairs.tsv
                 └─────────────────────────────────────────────┘──► synthetic_features.tsv
                 ┌─────────────────────────────────────────────┐
                 │ 2. 10 stratified train/test splits (70/30)  │
                 │    on each training part:                   │
                 │      MDFS 10 runs → consensus pairs         │
                 │      Baseline RF: original features         │──► auc_per_split.tsv
                 │      MDFS RF: original + LR_/GM_ features   │──► split_consensus_pairs.tsv
                 │    both scored on the test part             │──► selected_features.tsv
                 └─────────────────────────────────────────────┘
                 ┌─────────────────────────────────────────────┐
                 │ 3. ROC curves and summary                   │──► roc_comparison.png/.pdf
                 └─────────────────────────────────────────────┘──► summary.txt
```

### 1. MDFS pair detection (one run, `run_mdfs.R`)

1. **1D MDFS**: information gain (IG) of each feature about the label.
2. **2D MDFS**: features that carry significant information jointly with at
   least one partner.
3. **Pair enumeration** (`ComputeInterestingTuples`): for each pair
   (base, contributing) the IG the contributing feature adds on top of the
   base feature.
4. **Bonferroni-corrected pair-level threshold**: the added IG must exceed a
   chi-squared cutoff corrected for all n × (n − 1) / 2 feature pairs.
5. **True-synergy filter**: the pair must carry more information than either
   feature alone, `total_IG > max(base_IG_1D, contributing_IG_1D)`.

For a pair of features:

| Quantity | Meaning |
|---|---|
| `base_IG_1D` | IG of the base feature alone, I(X_base; Y) |
| `IG_2D_added` | IG the contributing feature adds given the base, I(X_contributing; Y \| X_base) |
| `total_IG` | joint IG of the pair, `base_IG_1D + IG_2D_added` |
| `ig_gain_pct` | `IG_2D_added / base_IG_1D × 100` |

Each pair is reported in both directions (base and contributing swapped).
IG values are as reported by the MDFS package; they are comparable within a
dataset, not across datasets of different size. `ig_gain_pct` becomes very
large when `base_IG_1D` is close to zero, so the table is sorted by
`total_IG` instead.

### 2. Consensus pairs (10 runs)

MDFS discretises features with a random component, so results vary between
runs. MDFS is therefore run **10 times** (seeds 42–51). In each run, unordered
pairs are scored by the larger `total_IG` of their two directions and the
**top 2 %** are selected (at least one pair). A pair selected in **at least 6
of 10 runs** is a **consensus pair**.

### 3. Synthetic features

For each consensus pair (f1, f2):

- `LR_f1__f2`: log-ratio, `log((f1 + ε) / (f2 + ε))`
- `GM_f1__f2`: geometric mean, `sqrt((f1 + ε) × (f2 + ε))`

with ε = 1e-9.

### 4. Model comparison

The data are split **10 times** into stratified 70 % training / 30 % test
parts. In each split, MDFS consensus pairs are found **on the training part
only** (10 MDFS runs, ≥ 6 occurrences), so the test samples never influence
pair selection. Two Random Forests are then trained on the same training part:

- **Baseline RF**: the original features
- **MDFS RF**: the original features plus the `LR_` / `GM_` features of the
  consensus pairs

Both use the same procedure. The number of features k (50, 100, …, 500, or
all features if there are fewer than 50) is tuned on an internal 75/25 split
of the training part. The top-k features by RF importance are selected, and a
final RF (1000 trees, `max_features = sqrt`, balanced class weights) is scored
on the test part. If a split has no consensus pairs, both models use the same
features and give the same result.

This is the procedure used in the manuscript's leave-one-cohort-out analysis
(10 MDFS runs, top 2 %, ≥ 6/10 consensus, same RF), applied to train/test
splits of a single cohort.

---

## Requirements

| Dependency | Version | Tested with |
|---|---|---|
| R | >= 4.0 | 4.5.1 |
| R package `MDFS` | >= 1.5 | 1.5.5 |
| R package `data.table` | any recent | 1.18.4 |
| Python | >= 3.8 | 3.9.6 |
| `pandas` | >= 1.3 | 2.3.1 |
| `numpy` | >= 1.21 | 2.0.2 |
| `scikit-learn` | >= 1.0 | 1.6.1 |
| `matplotlib` | >= 3.5 | 3.9.4 |

`Rscript` must be on your `PATH`.

```r
install.packages(c("MDFS", "data.table"))
```

```bash
pip install -r requirements.txt
```

---

## Usage

```bash
python run_analysis.py [--X FILE] [--y FILE] [--out-dir DIR] [options]
```

| Option | Default | Description |
|---|---|---|
| `--X` | `test_files/X_test.tsv` | Feature matrix |
| `--y` | `test_files/y_test.tsv` | Labels |
| `--out-dir` | `results` | Output directory |
| `--n-mdfs-runs` | 10 | MDFS runs per pair selection |
| `--min-runs` | 6 | Minimum runs in which a pair must be selected to be a consensus pair |
| `--top-fraction` | 0.02 | Fraction of pairs selected per MDFS run (at least one) |
| `--n-splits` | 10 | Stratified train/test splits |
| `--test-size` | 0.3 | Test fraction of each split |
| `--seed` | 42 | Seed of the splits; MDFS run *i* uses seed + *i* |

**Runtime.** MDFS runs 10 × (1 + number of splits) times, so 110 times with the
defaults. Per split and model, up to 13 Random Forests of 1000 trees are
trained (11 to tune k, 2 for the final model; 4 in total when there are fewer
than 50 features). The test data take about a minute. For a cohort with a few
hundred samples and a few thousand features, expect tens of minutes to a few
hours.

With few features there are few pairs, so the top 2 % may be a single pair
(ties at the cutoff are broken by feature name). Raise `--top-fraction` if no
consensus pairs are found.

### Input format

**Feature matrix (`--X`)**: tab-separated, samples in rows, first column is
the sample ID. Values must be numeric and non-negative, typically (relative)
abundances; log-ratio and geometric-mean features cannot be computed from
CLR- or log-transformed data. Missing values are not allowed. Spaces in feature
names are replaced by underscores.

```
sample_id   Fusobacterium_nucleatum   Parvimonas_micra   ...
S001        0.02831045                0.00112034         ...
```

**Labels (`--y`)**: tab-separated, with a header, two columns:

```
sample_id   label
S001        0
S151        1
```

`0` = disease / case, `1` = control / healthy. Only samples present in both
files are used; the number of samples found in only one file is reported.
Sample IDs must be unique. The runner stops with an explanatory message when
an input does not follow this format (e.g. a comma-separated file, a labels
file without header, labels other than 0/1).

### Outputs

| File | Content |
|---|---|
| `mdfs_synergies.tsv` | Full dataset: every (base, contributing) direction that passed both MDFS filters in at least one run, averaged over those runs. Columns: `base`, `contributing`, `base_IG_1D_mean`, `contributing_IG_1D_mean`, `IG_2D_added_mean`, `IG_2D_added_sd`, `total_IG_mean`, `ig_gain_pct`, `n_runs_significant` (runs in which the direction passed the filters), `n_runs_top` (runs in which the pair was in the top 2 %), `consensus` (`n_runs_top` ≥ 6). Sorted by `total_IG_mean`. |
| `consensus_pairs.tsv` | Full dataset: consensus pairs (`feature_1`, `feature_2`, `n_runs_top`, `total_IG_mean`). Here `total_IG_mean` is the pair's score (larger `total_IG` of its two directions) averaged over the runs in which it was in the top 2 %. |
| `synthetic_features.tsv` | Full dataset: original features plus `LR_` / `GM_` features of the consensus pairs, ready for your own models |
| `auc_per_split.tsv` | Per split and model: `n_train`, `n_test`, `n_consensus_pairs`, `n_features`, `best_k`, `auc`, `accuracy` |
| `split_consensus_pairs.tsv` | `feature_1`, `feature_2`, `n_splits`: number of splits in which each pair was a consensus pair |
| `selected_features.tsv` | `feature`, `n_splits_selected`, `synthetic`: MDFS RF, number of splits in which each feature was among the top-k |
| `roc_comparison.png` / `.pdf` | Left: mean ROC curve ± sd of both models over the test splits, with individual splits faint. Right: test AUC of both models in each split. |
| `summary.txt` | Settings, number of pairs, mean ± sd AUC and accuracy of both models, and the per-split AUC difference |

The AUCs come from the held-out test parts. The consensus pairs in
`mdfs_synergies.tsv` / `consensus_pairs.tsv` come from the full dataset, so
use them to describe the cohort, not to estimate performance. The 10 test
parts overlap, so the per-split AUCs are not independent. AUC and ROC curves
use label 1 (control) as the positive class; AUC is the same either way.

---

## Test data

`test_files/X_test.tsv` and `test_files/y_test.tsv` are a synthetic
CRC-like dataset (300 samples: 150 disease, 150 control; 30 species as
relative abundances), generated by `test_files/make_test_data.py`:

- **Synergistic pair**: *Parvimonas micra* and *Gemella morbillorum* share a
  strongly varying sample-specific abundance, so each is only weakly
  informative alone, but their ratio is high in disease and low in controls.
- **Weak individual markers**: *Fusobacterium nucleatum* and
  *Peptostreptococcus anaerobius* slightly higher in disease,
  *Faecalibacterium prausnitzii* and *Roseburia intestinalis* slightly lower.
- The other 24 species are noise.

```bash
python run_analysis.py
```

Expected result (`test_files/expected/`):

- `Gemella_morbillorum` / `Parvimonas_micra` has the highest `total_IG` and
  is the only consensus pair, both on the full dataset and in all 10 splits.
- Its `LR_` and `GM_` features are included in the MDFS RF in every split
  (with 30 features, k-tuning keeps all features, so this is expected).
- Test AUC: Baseline RF 0.863 ± 0.025, MDFS RF 0.933 ± 0.026; MDFS RF is
  higher in 10/10 splits.

Exact AUCs can differ slightly between scikit-learn versions.

---

## Repository structure

```
├── run_analysis.py          # Runner: the whole analysis, outputs and figure
├── mdfs_pairs.py            # Runs run_mdfs.R repeatedly; synergy table and consensus pairs
├── rf_model.py              # LR_/GM_ synthetic features; Random Forest with top-k selection
├── run_mdfs.R               # One MDFS run: 1D/2D MDFS, Bonferroni threshold, synergy filter
├── requirements.txt
└── test_files/
    ├── X_test.tsv, y_test.tsv   # Synthetic test data
    ├── make_test_data.py        # Generates the test data
    └── expected/                # summary.txt and roc_comparison.png of the test run
```

`run_mdfs.R` must stay in the same directory as the Python scripts.
