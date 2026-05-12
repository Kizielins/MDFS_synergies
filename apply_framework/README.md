# Apply MDFS framework to other cohorts

Apply the MDFS framework to new cohorts: run 1D/2D MDFS (3 runs per cohort, aggregate with mean) and report information gains; load cohort data (taxa + functions combined) and build tables of feature pairs showing 2D IG gain over the base feature's 1D IG.

## Labels

- **1 (control/healthy):** `study_condition` or `disease` = `control` or `healthy` (case-insensitive).
- **0 (disease):** Any other value (CRC, Adenoma, or any other disease).
Samples with missing/unrecognized labels are dropped.

## Data layout

Put input data in `data/`. Two layouts are supported:

1. **Merged**
   - `data/merged_metadata.tsv` — columns: `sample_id`, `study_name` or `source_cohort`, and `study_condition` or `disease`.
   - `data/merged_abundances.txt` — taxa (rows = species, columns = sample IDs).
   - `data/merged_pathways.txt` — pathways (rows = pathway names, columns = sample IDs).
   Cohorts = unique `study_name` or `source_cohort`.

2. **Per-cohort**
   - `data/metadata/<cohort>.tsv`, `data/taxa/<cohort>.tsv`, `data/pathways/<cohort>.tsv`.

## MDFS 2D IG interpretation

The 2D information gain from MDFS is **additional** IG. For a pair (base, contributing):
- `base_IG_1D`: the base feature's 1D information gain
- `IG_2D_added`: how much the contributing variable improves the base feature's IG
- `total_IG = base_IG_1D + IG_2D_added`

Each pair is reported in **both directions** (base/contributing are asymmetric).

## Scripts

### 1. `run_mdfs_1d_2d.py` — run MDFS 1D and 2D, report IGs for pairs

Runs the repo's `run_mdfs.R` **3 times per cohort** (like the main pipeline), then aggregates IGs with **mean** across runs and writes:

- `mdfs_2d_<cohort>.tsv` — 2D tuples: `base`, `contributing`, `base_IG_1D`, `contributing_IG_1D`, `IG_2D_added`, `total_IG` (from last run)
- `mdfs_1d_ig_<cohort>.tsv` — per-feature 1D IG: `Feature`, `IG`, `p_value`, `adjusted_p_value` (from last run)
- `mdfs_ig_info_<cohort>.tsv` — one row per **run** per (base, contributing) pair: `cohort`, `run`, `base`, `contributing`, `base_IG_1D`, `contributing_IG_1D`, `IG_2D_added`, `total_IG`
- `ig_gain_summary_<cohort>.tsv` — one row per pair, **mean across 3 runs**: `base_IG_1D_mean`, `contributing_IG_1D_mean`, `IG_2D_added_mean`, `IG_2D_added_std`, `total_IG_mean`, `n_runs`, `ig_gain`, `ig_gain_pct`

**Examples:**

```bash
cd apply_framework

# One cohort (from apply_framework/data)
python run_mdfs_1d_2d.py --data-dir data --cohort VogtmannE_2016 --out-dir results

# All cohorts in data (process all cohorts)
python run_mdfs_1d_2d.py --data-dir data --all-cohorts --out-dir results

# Pre-built X and y TSVs
python run_mdfs_1d_2d.py --X path/to/X.tsv --y path/to/y.tsv --cohort-name my_cohort --out-dir results
```

From repo root:

```bash
python apply_framework/run_mdfs_1d_2d.py --data-dir apply_framework/data --all-cohorts --out-dir apply_framework/results
```

Requires R and the `MDFS` package; `run_mdfs.R` must be in the repository root.

### 2. `load_cohort_data.py` — load cohorts, combined matrix, synergy table

- Discovers cohorts from `data/` (merged metadata or `metadata/*.tsv`).
- For each cohort: loads taxa + pathways, builds one combined matrix (samples x features), and gets disease labels from metadata.
- Reads `mdfs_ig_info_*.tsv` from the results dir and builds a **synergy table**: aggregated by (cohort, base, contributing) with mean across runs, computing `ig_gain` (= `IG_2D_added_mean`) and `ig_gain_pct` (= `ig_gain / base_IG_1D_mean * 100`).

**Examples:**

```bash
cd apply_framework

# List cohorts in data
python load_cohort_data.py --data-dir data --list-cohorts

# Build synergy table from existing MDFS results
python load_cohort_data.py --data-dir data --mdfs-results-dir results --out-table synergy_by_cohort.tsv
```

Output table columns: `cohort`, `base`, `contributing`, `base_IG_1D_mean`, `contributing_IG_1D_mean`, `IG_2D_added_mean`, `total_IG_mean`, `n_runs`, `ig_gain`, `ig_gain_pct`.

## How to process all cohorts

1. Put cohort data in `apply_framework/data/` (merged or per-cohort layout).
2. From the `apply_framework` directory, run MDFS for **all** cohorts (3 runs each, IGs aggregated with mean):

   ```bash
   cd apply_framework
   python run_mdfs_1d_2d.py --data-dir data --all-cohorts --out-dir results
   ```

   Or from the repository root:

   ```bash
   python apply_framework/run_mdfs_1d_2d.py --data-dir apply_framework/data --all-cohorts --out-dir apply_framework/results
   ```

3. Build the synergy table (one row per cohort per (base, contributing) pair, aggregated across runs):

   ```bash
   cd apply_framework
   python load_cohort_data.py --data-dir data --mdfs-results-dir results --out-table synergy_by_cohort.tsv
   ```
