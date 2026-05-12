# How to Run TAXA Analysis

## Prerequisites
1. Make sure `prep_data.ipynb` has been run to prepare the data files in `prepared_data/` directory
2. Ensure `framework.ipynb` has been executed (or it will auto-load when running `main_code.ipynb`)
3. Make sure R and the MDFS package are installed

## How to Run

1. **Open `main_code.ipynb` in Jupyter**
2. **Run all cells** (or use "Run All" from the Cell menu)
   - The notebook will automatically load functions from `framework.ipynb`
   - It will run the LOCO analysis for **TAXA only**

## Expected Outputs

The analysis will generate the following files in `results_taxa/`:

### 1. Performance Results
- **`loco_results_taxa.csv`**: AUC scores for each test cohort
  - Columns: Test Cohort, N Samples, All Features, Synthetic Only, All Features + Synthetic
  - Includes Mean and Weighted Mean rows
  
- **`loco_accuracy_taxa.csv`**: Accuracy scores for each test cohort
  - Same structure as AUC file but with accuracy values

### 2. Feature Selection Results
- **`base_features_taxa.csv`**: Base (non-synthetic) features selected by RF
  - Columns: Feature, Selection Count, Selection %
  - Shows how many folds each feature appeared in
  
- **`synthetic_features_taxa.csv`**: Synthetic features (LR_, GM_) selected by RF
  - Columns: Feature, Selection Count, Selection %
  - Shows synthetic features that were consistently selected

### 3. NEW: Consensus-Based Results (from 3 runs per fold)
- **`mdfs_ig_info_taxa.csv`**: MDFS Information Gain values for all runs
  - Columns: fold, test_cohort, run, feature1, feature2, pair_ig, f1_ig, f2_ig
  - Contains IG values for:
    - Pair IG: 2D information gain for feature pairs
    - f1_ig, f2_ig: 1D information gain for individual features
  - **Purpose**: Allows comparison to show pair_ig > individual_ig (synergy)
  
- **`rf_features_per_fold_taxa.csv`**: RF selected features for each run
  - Columns: fold, test_cohort, run, n_features, features
  - Shows which features were selected in each of the 3 runs per fold
  - **Purpose**: Track feature stability across runs

## What Happens During Execution

1. **For each LOCO fold** (14 test cohorts):
   - Runs **3 times** with different random seeds (42, 43, 44)
   - Each run:
     - Runs MDFS to find top 20% feature pairs (by IG)
     - Generates synthetic features (LR_, GM_)
     - Trains RF and selects features
     - Saves IG values and selected features
   
2. **Consensus selection**:
   - MDFS pairs: Keeps pairs appearing in >= 2 out of 3 runs
   - RF features: Tracks features appearing in >= 2 out of 3 runs
   
3. **Final evaluation**:
   - Uses consensus pairs to generate synthetic features
   - Trains final RF model and evaluates on test set
   - Reports AUC and accuracy

## Console Output

You'll see:
- Progress for each fold (1/14, 2/14, etc.)
- For each fold: 3 runs (RUN 1/3, RUN 2/3, RUN 3/3)
- Consensus information (how many pairs/features in consensus)
- Performance metrics (AUC, accuracy) for each scenario
- Summary statistics at the end

## Expected Runtime

- **~3x longer** than before (since each fold runs 3 times)
- For 14 cohorts with 3 runs each = 42 total runs
- Each run includes: MDFS computation + RF training
- Estimated: Several hours depending on data size and hardware
