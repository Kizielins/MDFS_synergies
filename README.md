# Bacterial Feature Synergy Analysis

Identifies synergistic pairs of bacterial (or functional) features using
**Multidimensional Feature Selection (MDFS)** information gain.

For every feature pair returned by MDFS the tool reports:
- 1D information gain (IG) for each feature individually
- 2D information gain for the pair jointly
- **% synergy gain** = `(pair_IG − median(IG_f1, IG_f2)) / median(IG_f1, IG_f2) × 100`

MDFS is run multiple times with different random seeds; IGs are averaged across
runs for robustness.  The output table is sorted by % synergy gain (descending).

---

## Requirements

| Dependency | Version |
|---|---|
| R | ≥ 4.0 |
| R package `MDFS` | ≥ 1.5 |
| R package `data.table` | any recent |
| Python | ≥ 3.8 |
| Python package `pandas` | any recent |

### Install R packages

```r
install.packages("data.table")
install.packages("MDFS")
```

---

## Files

```
bacterial_synergy_analysis/
├── run_mdfs.R              # R script: computes 1D and 2D MDFS information gains
├── compute_synergies.py    # Python script: orchestrates runs, computes synergy metrics, writes output
└── test_files/
    ├── X_test.tsv          # 40-sample × 12-feature synthetic microbiome matrix
    └── y_test.tsv          # Labels (0 = disease, 1 = control)
```

Both scripts must reside in the **same directory** — `compute_synergies.py` locates
`run_mdfs.R` relative to itself.

---

## Usage

```
python compute_synergies.py --X <features.tsv> --y <labels.tsv> [--out results.txt] [--n-runs N]
```

| Argument | Description | Default |
|---|---|---|
| `--X` | Feature matrix TSV (samples × features, first col = sample IDs) | required |
| `--y` | Labels TSV (two cols: sample_id, label; `0`=disease, `1`=control/healthy) | required |
| `--out` | Output file path | stdout |
| `--n-runs` | Number of MDFS runs to average (uses seeds 12–16) | 3 |

### Quick start

```bash
python compute_synergies.py \
    --X test_files/X_test.tsv \
    --y test_files/y_test.tsv \
    --out test_files/results.txt
```

Progress is printed to stderr; the table goes to stdout (or `--out`).

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

Tab-separated plain-text table, sorted by % synergy gain (highest first):

```
#   Feature_1               Feature_2           IG_F1       IG_F2       IG_Pair     Pct_Synergy_Gain
1   Parvimonas_micra        Gemella_morbillorum  0.052341    0.049817    0.381204    668.12
2   Fusobacterium_nucleatum Peptostrep_anaerobius 0.412300   0.183600    0.521400    65.30
...
```

| Column | Description |
|---|---|
| `#` | Row number |
| `Feature_1` | First feature in the pair |
| `Feature_2` | Second feature in the pair |
| `IG_F1` | Mean 1D information gain for Feature_1 (averaged across runs) |
| `IG_F2` | Mean 1D information gain for Feature_2 (averaged across runs) |
| `IG_Pair` | Mean 2D information gain for the pair (averaged across runs) |
| `Pct_Synergy_Gain` | `(IG_Pair − median(IG_F1, IG_F2)) / median(IG_F1, IG_F2) × 100` |

A positive `Pct_Synergy_Gain` means the pair carries **more** discriminative
information together than either feature does alone.

---

## Test data

`test_files/X_test.tsv` and `test_files/y_test.tsv` contain a synthetic CRC-like
microbiome dataset (40 samples, 12 bacterial features).  The data is designed so
that **Parvimonas_micra** and **Gemella_morbillorum** form a synergistic pair:
in disease samples, one or the other is elevated (never both), while in controls
both are low — so neither bacterium alone is strongly discriminative but the pair
jointly is.

Run the test:

```bash
python compute_synergies.py \
    --X test_files/X_test.tsv \
    --y test_files/y_test.tsv \
    --out test_files/results.txt

cat test_files/results.txt
```

Expected: `Parvimonas_micra / Gemella_morbillorum` appears near the top with the
highest % synergy gain.

---

## How it works

1. **MDFS 2D** (`ComputeMaxInfoGains(..., dimensions=2)`) — identifies all feature
   pairs and their joint information gain with the class label.
2. **MDFS 1D** (`ComputeMaxInfoGains(..., dimensions=1)`) — computes per-feature
   individual information gain.
3. Steps 1–2 are repeated `--n-runs` times with different random seeds; IGs are
   averaged across runs.
4. For each pair the **% synergy gain** is computed and the table is sorted
   accordingly.
