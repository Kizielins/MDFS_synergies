# Metagenomic Feature Synergy Analysis

Identifies synergistic pairs of metagenomic features using
**Multidimensional Feature Selection (MDFS)** information gain in 1D and 2D,
and turns them into synthetic features for classification.

The tool works on a **single cohort**: one feature matrix and one set of
binary labels.

| Script | What it does |
|---|---|
| `compute_synergies.py` | Finds synergistic feature pairs and reports their information gains |
| `synthetic_features.py` | Builds log-ratio / geometric-mean features from the top synergistic pairs and evaluates them with a Random Forest (cross-validation) |
| `run_mdfs.R` | MDFS engine called by both scripts (not run directly) |

## How MDFS 2D information gain works

MDFS evaluates each feature's individual discriminative power (1D) and then
tests whether pairing two features yields additional predictive information
(2D). For a pair of features designated **base** and **contributing**:

- **`base_IG_1D`** — the base feature's individual information gain about the
  class label, i.e. I(X_base; Y).
- **`IG_2D_added`** — the *additional* IG that the contributing feature brings
  when combined with the base feature, i.e. I(X_contributing; Y | X_base).
  This is the conditional mutual information of the contributing feature given
  the base feature.
- **`total_IG`** = `base_IG_1D + IG_2D_added` — the joint information gain of
  the pair, i.e. I(X_base, X_contributing; Y).

Each pair is reported in **both directions** (base/contributing roles are
asymmetric: the decomposition changes, but `total_IG` remains the same).

The **% synergy gain** (`ig_gain_pct`) measures how much the contributing
feature amplifies the base feature's signal:

```
ig_gain_pct = IG_2D_added / base_IG_1D × 100
```

A value of 100% means the pair's joint IG is double the base feature's
individual IG; values above 100% indicate that the contributing feature adds
more information than the base carried alone.

### Statistical thresholds

Two filters ensure that reported pairs represent genuine synergies:

1. **Bonferroni-corrected pair-level threshold** — `ComputeInterestingTuples`
   in R uses an IG cutoff derived from a chi-squared quantile corrected for the
   total number of feature pairs tested (n × (n − 1) / 2 for n features).
   Only pairs exceeding this threshold are retained.

2. **True-synergy filter** — the pair's `total_IG` must exceed the larger of
   the two individual 1D IGs (`total_IG > max(base_IG_1D, contributing_IG_1D)`).
   This removes artifacts where a strong feature is paired with noise.

MDFS is run **3 times** with different random seeds (`--n-runs`); IGs are
averaged across runs for robustness.

---

## Requirements

| Dependency | Version | Tested with |
|---|---|---|
| R | >= 4.0 | 4.5.1 |
| R package `MDFS` | >= 1.5 | 1.5.5 |
| R package `data.table` | any recent | 1.18.4 |
| Python | >= 3.8 | 3.10.13 |
| Python package `pandas` | >= 1.3 | 2.3.2 |
| Python package `numpy` | >= 1.21 | 1.26.4 |
| Python package `scikit-learn` | >= 1.0 (only for `synthetic_features.py`) | 1.7.2 |

`Rscript` must be on your `PATH`.

### Install

R packages:

```r
install.packages(c("MDFS", "data.table"))
```

Python packages:

```bash
pip install -r requirements.txt
```

---

## Repository structure

```
MDFS_synergies/
├── compute_synergies.py     # Synergistic pairs and their information gains
├── synthetic_features.py    # Synthetic features from top pairs + Random Forest evaluation
├── run_mdfs.R               # MDFS 1D/2D, Bonferroni threshold, true-synergy filter
├── requirements.txt         # Python dependencies
└── test_files/
    ├── X_test.tsv           # 40-sample × 12-feature synthetic dataset
    ├── y_test.tsv           # Labels (0 = disease, 1 = control)
    └── expected/            # Expected outputs of the test commands below
        ├── synergies.tsv
        └── synthetic/
```

Run all commands from the repository root. `compute_synergies.py` and
`synthetic_features.py` must stay in the same directory as `run_mdfs.R`.

---

## Usage

### 1. Synergistic pairs — `compute_synergies.py`

```bash
python compute_synergies.py --X my_features.tsv --y my_labels.tsv --out synergies.tsv
```

| Option | Default | Description |
|---|---|---|
| `--X` | required | Feature matrix (see [Input format](#input-format)) |
| `--y` | required | Labels |
| `--out` | stdout | Output file |
| `--n-runs` | 3 | MDFS runs to average over (max 5) |

Output: one row per (base, contributing) pair, sorted by `ig_gain_pct`
(descending):

| Column | Description |
|---|---|
| `#` | Rank |
| `base` | Base feature in the pair |
| `contributing` | Contributing feature in the pair |
| `base_IG_1D` | Mean 1D IG of the base feature across runs |
| `contributing_IG_1D` | Mean 1D IG of the contributing feature across runs |
| `IG_2D_added` | Mean additional IG from the contributing feature (conditional on base) |
| `total_IG` | Mean joint IG of the pair (`base_IG_1D + IG_2D_added`) |
| `ig_gain_pct` | `IG_2D_added / base_IG_1D × 100` |

A pair is averaged over the runs in which it passed both filters.

### 2. Synthetic features — `synthetic_features.py`

```bash
python synthetic_features.py --X my_features.tsv --y my_labels.tsv --out-dir out
```

For each selected pair (f1, f2) two synthetic features are created:

- `LR_f1__f2` — log-ratio, `log((f1 + ε) / (f2 + ε))`
- `GM_f1__f2` — geometric mean, `sqrt((f1 + ε) × (f2 + ε))`

**Pair selection.** MDFS is run `--n-runs` times. In each run, unordered pairs
are scored by `total_IG` and the top `--top-fraction` are kept (at least one).
Pairs kept in **every** run form the consensus set.

**Evaluation.** Stratified k-fold cross-validation (`--n-folds`). Pair
selection is repeated inside each training fold, so the test fold never
informs which pairs are used. Three feature sets are compared:

- `All Features` — original features
- `Synthetic Only` — `LR_` / `GM_` features
- `All Features + Synthetic` — both

For each feature set, the number of features *k* (50–500, or all if fewer) is
tuned on an internal 75/25 split of the training fold, the top-*k* features by
Random Forest importance are selected, and a final Random Forest (1000 trees,
balanced class weights) is scored on the test fold. Feature sets with fewer
than 10 features are not evaluated (reported as empty).

| Option | Default | Description |
|---|---|---|
| `--X` | required | Feature matrix |
| `--y` | required | Labels |
| `--out-dir` | required | Output directory |
| `--n-runs` | 3 | MDFS runs per pair selection (max 5) |
| `--top-fraction` | 0.01 | Fraction of pairs kept per MDFS run |
| `--n-folds` | 5 | Cross-validation folds |
| `--no-cv` | off | Skip evaluation; only select pairs and write synthetic features |

Outputs in `--out-dir`:

| File | Content |
|---|---|
| `synergy_pairs.tsv` | Consensus pairs selected on the full dataset: `feature_1`, `feature_2`, `total_IG_mean` |
| `synthetic_features.tsv` | Original features plus `LR_` / `GM_` features for those pairs (samples × features), ready for your own models |
| `cv_performance.tsv` | Per fold and feature set: `n_test`, `n_pairs`, `n_features`, `best_k`, `auc`, `accuracy` |
| `selected_features.tsv` | For `All Features + Synthetic`: how many folds selected each feature, and whether it is synthetic |

A summary (mean ± sd AUC and accuracy per feature set) is printed at the end.

The default `--top-fraction 0.01` suits datasets with hundreds or thousands of
features. With few features (few pairs) raise it, otherwise the consensus set
may be empty.

---

## Input format

**Feature matrix (`--X`)**

Tab-separated, samples in rows, first column is the sample ID. Values are
typically relative abundances (taxa, pathways, or both combined):

```
sample_id   Fusobacterium_nucleatum   Peptostreptococcus_anaerobius   ...
S001        0.28310452                0.11203471                      ...
S002        0.01823940                0.03941200                      ...
```

**Labels (`--y`)**

Tab-separated, two columns with a header row:

```
sample_id   label
S001        0
S002        0
S021        1
```

`0` = disease / case, `1` = control / healthy. Only samples present in both
files are used. Spaces in feature names are replaced by underscores.

---

## Test data

`test_files/X_test.tsv` and `test_files/y_test.tsv` contain a synthetic CRC-like
microbiome dataset (40 samples, 12 bacterial features). The data is designed so
that **Parvimonas_micra** and **Gemella_morbillorum** form a synergistic pair:
in disease samples, one or the other is elevated (never both), while in controls
both are low — so neither bacterium alone is strongly discriminative but the pair
jointly is.

```bash
python compute_synergies.py --X test_files/X_test.tsv --y test_files/y_test.tsv --out out/synergies.tsv
python synthetic_features.py --X test_files/X_test.tsv --y test_files/y_test.tsv --out-dir out/synthetic --top-fraction 0.2
```

Compare with the expected results:

```bash
diff out/synergies.tsv test_files/expected/synergies.tsv
diff -r out/synthetic test_files/expected/synthetic
```

Expected: both `Parvimonas_micra / Gemella_morbillorum` and
`Gemella_morbillorum / Parvimonas_micra` pass both synergy filters, with
`ig_gain_pct` ≈ 1076% and ≈ 396%. The small dataset is a check that the tools
run: several features separate the classes on their own, so the
cross-validated AUC is 1.0 for every evaluated feature set and
`Synthetic Only` has too few features to be evaluated.
