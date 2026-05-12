# Metagenomic Feature Synergy Analysis

Identifies synergistic pairs of metagenomic features using
**Multidimensional Feature Selection (MDFS)** information gain in 1D and 2D.

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
   total number of feature pairs tested (~1.5 million for 1,738 features).
   Only pairs exceeding this threshold are retained.

2. **True-synergy filter** — the pair's `total_IG` must exceed the larger of
   the two individual 1D IGs (`total_IG > max(base_IG_1D, contributing_IG_1D)`).
   This removes artifacts where a strong feature is paired with noise.

MDFS is run **3 times** with different random seeds per cohort; IGs are
averaged across runs for robustness.

---

## Requirements

| Dependency | Version |
|---|---|
| R | >= 4.0 |
| R package `MDFS` | >= 1.5 |
| R package `data.table` | any recent |
| Python | >= 3.8 |
| Python package `pandas` | any recent |

### Install R packages

```r
install.packages("data.table")
install.packages("MDFS")
```

---

## Repository structure

```
MDFS_synergies/
├── run_mdfs.R                  # R script: 1D and 2D MDFS, Bonferroni threshold,
│                               #   synergy filter, ComputeInterestingTuples
├── compute_synergies.py        # Python: orchestrates runs, computes synergy metrics
├── run_taxa_analysis.py        # Python: CRC taxa-only analysis pipeline
├── test_files/
│   ├── X_test.tsv              # 40-sample × 12-feature synthetic dataset
│   └── y_test.tsv              # Labels (0 = disease, 1 = control)
├── apply_framework/            # Apply the framework to non-CRC cohorts
│   ├── run_mdfs_1d_2d.py       # Run MDFS 1D+2D per cohort (3 runs, mean aggregation)
│   ├── load_cohort_data.py     # Load cohort data, build synergy table
│   ├── analyse_synergy_gains.ipynb  # Analysis notebook: figures, tables, markers
│   └── results/                # Per-cohort MDFS outputs
├── results_combined/           # CRC pipeline results (taxa + functions combined)
├── results_functions/          # CRC pipeline results (functions only)
└── results_taxa/               # CRC pipeline results (taxa only)
```

---

## Usage

### Quick start (test data)

```bash
python compute_synergies.py \
    --X test_files/X_test.tsv \
    --y test_files/y_test.tsv \
    --out test_files/results.txt
```

### Apply framework to new cohorts

See [`apply_framework/README.md`](apply_framework/README.md) for full details.

```bash
cd apply_framework

# Run MDFS for all cohorts
python run_mdfs_1d_2d.py --data-dir data --all-cohorts --out-dir results

# Build synergy table
python load_cohort_data.py --data-dir data --mdfs-results-dir results --out-table synergy_by_cohort.tsv
```

---

## Input format

**Feature matrix (`--X`)**

Tab-separated, first column is the sample ID:

```
sample_id   Fusobacterium_nucleatum   Peptostreptococcus_anaerobius   ...
S001        0.28310452                0.11203471                      ...
S002        0.01823940                0.03941200                      ...
```

**Labels (`--y`)**

Tab-separated, two columns:

```
sample_id   label
S001        0
S002        0
S021        1
```

`0` = disease / case, `1` = control / healthy.

---

## Output format

The per-cohort output (`ig_gain_summary_<cohort>.tsv`) is a tab-separated table
with one row per (base, contributing) pair, averaged across 3 MDFS runs:

| Column | Description |
|---|---|
| `base` | Base feature in the pair |
| `contributing` | Contributing feature in the pair |
| `base_IG_1D_mean` | Mean 1D IG of the base feature across runs |
| `contributing_IG_1D_mean` | Mean 1D IG of the contributing feature across runs |
| `IG_2D_added_mean` | Mean additional IG from the contributing feature (conditional on base) |
| `IG_2D_added_std` | Std dev of `IG_2D_added` across runs |
| `total_IG_mean` | Mean joint IG of the pair (`base_IG_1D + IG_2D_added`) |
| `n_runs` | Number of runs in which this pair was significant |
| `ig_gain` | Same as `IG_2D_added_mean` |
| `ig_gain_pct` | `IG_2D_added_mean / base_IG_1D_mean × 100` |

---

## Test data

`test_files/X_test.tsv` and `test_files/y_test.tsv` contain a synthetic CRC-like
microbiome dataset (40 samples, 12 bacterial features). The data is designed so
that **Parvimonas_micra** and **Gemella_morbillorum** form a synergistic pair:
in disease samples, one or the other is elevated (never both), while in controls
both are low — so neither bacterium alone is strongly discriminative but the pair
jointly is.

Expected: `Parvimonas_micra / Gemella_morbillorum` appears near the top with the
highest % synergy gain.
