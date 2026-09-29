# Metagenomic Feature Synergy Analysis

Finds synergistic pairs of metagenomic features (taxa, pathways, or both) with
**Multidimensional Feature Selection (MDFS)**, turns the most robust pairs into
synthetic features, and tests whether they improve a Random Forest classifier.

The tool analyses a **single cohort**: one feature matrix and one set of
binary labels (disease vs control).

```bash
python run_analysis.py                                   # bundled test data
python run_analysis.py --X my_features.tsv --y my_labels.tsv --out-dir my_results
python run_analysis.py --X my_features.tsv --y my_labels.tsv --no-evaluation   # synergies only, fast
```

> **Labels: `0` = disease / case, `1` = control / healthy.** This is the
> reverse of the common 1 = case convention; recode your labels if needed.

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

MDFS discretizes each feature at random. Within one run all steps use the
same seed and discretization, so the 1D and 2D IGs are directly comparable and
`total_IG` is the same in both directions of a pair.

For a pair of features:

| Quantity | Meaning |
|---|---|
| `base_IG_1D` | IG of the base feature alone, I(X_base; Y) |
| `IG_2D_added` | IG the contributing feature adds given the base, I(X_contributing; Y \| X_base) |
| `total_IG` | joint IG of the pair, `base_IG_1D + IG_2D_added` |
| `ig_gain_pct` | `IG_2D_added / base_IG_1D × 100` |

Each pair is reported in both directions (base and contributing swapped).
IG values are as reported by the MDFS package, which scales information gain
by the number of samples (IG in bits × n). This is why values are much larger
than 1 bit, and why they are comparable within a dataset but not across
datasets of different size. `ig_gain_pct` becomes very
large when `base_IG_1D` is close to zero, so the table is sorted by
`total_IG` instead.

### 2. Consensus pairs (10 runs)

Because of the random discretization, results vary between runs. MDFS is
therefore run **10 times** (seeds 42–51). In each run, unordered
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

---

## Requirements

| Dependency | Version | Tested with |
|---|---|---|
| R | >= 4.0 | 4.5.1 |
| R package `MDFS` | >= 1.5 | 1.5.5 |
| R package `data.table` | any recent | 1.18.4 |
| Python | >= 3.9 | 3.9.6, 3.10.13 |
| `pandas` | >= 1.4 | 1.4.0, 2.3.3 |
| `numpy` | >= 1.21 | 1.21.6, 2.2.6 |
| `scikit-learn` | >= 1.0 | 1.0.2, 1.7.2 |
| `matplotlib` | >= 3.5 | 3.5.3, 3.10.9 |

Tested on macOS (Apple silicon). `Rscript` must be on your `PATH`; the runner
checks for R, MDFS and data.table at startup and tells you what is missing.
MDFS's optional CUDA support is not used.

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
| `--seed` | 42 | Seed of the splits and Random Forests; MDFS run *i* uses seed + *i* |
| `--no-evaluation` | off | Only step 1 (synergies on the full dataset); skip the Random Forest comparison |

**Runtime.** MDFS runs 10 × (1 + number of splits) times, so 110 times with the
defaults (10 times with `--no-evaluation`). Per split and model, up to 13 Random Forests of 1000 trees are
trained (11 to tune k, 2 for the final model; 4 in total when there are fewer
than 50 features). The test data take about a minute. For a cohort with a few
hundred samples and a few thousand features, expect tens of minutes to a few
hours; use `--no-evaluation` if you only need the synergy table.

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

**Preprocessing.** The tool does not filter features. Before running it, we
recommend:

- removing rare features, e.g. those present in fewer than 10 % of samples:
  mostly-zero features give unstable IGs and many spurious pairs, and the
  number of pairs grows with the square of the number of features;
- using relative (or absolute) abundances, not CLR/log-transformed values;
- having at least a few dozen samples per class; with small cohorts, MDFS may
  find no significant pairs.

### Outputs

| File | Content |
|---|---|
| `mdfs_synergies.tsv` | Full dataset: every (base, contributing) direction that passed both MDFS filters in at least one run, averaged over those runs. Columns: `base`, `contributing`, `base_IG_1D_mean`, `contributing_IG_1D_mean`, `IG_2D_added_mean`, `IG_2D_added_sd`, `total_IG_mean`, `ig_gain_pct`, `n_runs_significant` (runs in which this direction passed the filters), `n_runs_top` (runs in which the unordered pair, in either direction, was in the top 2 %; can exceed `n_runs_significant`), `consensus` (`n_runs_top` ≥ 6). Sorted by `total_IG_mean`. Means are over the runs in which the direction passed, so the two directions of a pair can differ when their `n_runs_significant` differ. |
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

With `--no-evaluation` only `mdfs_synergies.tsv`, `consensus_pairs.tsv`,
`synthetic_features.tsv` and `summary.txt` are written. Existing files in the
output directory are overwritten (the runner prints a note), but files from an
earlier full run are not deleted, so use a new `--out-dir` to keep runs apart.

---

## Test data

`test_files/X_test.tsv` and `test_files/y_test.tsv` are a synthetic
CRC-like dataset (300 samples: 150 disease, 150 control; 30 species as
relative abundances), generated by `test_files/make_test_data.py`:

- **Synergistic pair**: *Parvimonas micra* and *Gemella morbillorum* share a
  strongly varying sample-specific abundance that blurs their individual
  signal, while their ratio is high in disease and low in controls. Each is
  somewhat informative alone; the pair is much more informative.
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
├── LICENSE
└── test_files/
    ├── X_test.tsv, y_test.tsv   # Synthetic test data
    ├── make_test_data.py        # Generates the test data
    └── expected/                # summary.txt and roc_comparison.png of the test run
```

`run_mdfs.R` must stay in the same directory as the Python scripts.

---

## Citation

If you use this tool, please cite:

- the accompanying manuscript (citation will be added on publication), and
- the MDFS package: Piliszek R, Mnich K, Migacz S, Tabaszewski P, Sułecki A,
  Polewko-Klim A, Rudnicki W. MDFS: MultiDimensional Feature Selection in R.
  *The R Journal* 11(1):198–210, 2019. doi:10.32614/RJ-2019-019

## License and contact

MIT License (see `LICENSE`). Questions and bug reports:
[GitHub issues](https://github.com/Kizielins/MDFS_synergies/issues).
