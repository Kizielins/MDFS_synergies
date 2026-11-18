# MDFS synergistic feature pairs in CRC and other microbiome cohorts — manuscript analysis code

This branch (`manuscript-revision`) contains the exact code and result tables behind every figure,
table and number of the manuscript *"Synergy-oriented feature selection …"* (Zielińska, Rudnicki, Łabaj;
npj Biofilms and Microbiomes, revision). The stand-alone, single-cohort version of the method is the
tool on the [`main`](https://github.com/Kizielins/MDFS_synergies) branch.

Folders are named after the manuscript sections they reproduce (`M` = Materials and Methods,
`R` = Results, in manuscript order). Each notebook starts with a header that says which figures,
tables and quoted numbers it produces, and every manuscript item is written to `manuscript/`.

## Where each item of the revised paper comes from

Figures and tables were renumbered late in the revision. The working file names in `manuscript/`
keep the earlier numbering, so use this table to go from the paper to the code. The files exactly as
submitted are in `manuscript/final/`.

| Paper | Submitted file (`manuscript/final/`) | Produced by | Working output |
|---|---|---|---|
| Figure 1 (a schematic, b all directed pairs) | `Figure_1.{pdf,png}` | `make_all_figures.py` ← `R1_synergistic_pairs/figure1_combined/make_figure1.py`, `figure1b_all_pairs/` | `R1_synergistic_pairs/figure1_combined/figure1_combined.png` |
| Figure 2 (hub-pair selection frequency, per-cohort effect sizes) | `Figure_2.{pdf,png}` | `make_all_figures.py` ← `R2_classification_performance/make_figure2.py` (`effect_sizes.ipynb`, `stable_pairs.ipynb`) | `R2_classification_performance/figure2.png` |
| Figure 3 (LR vs GM encodings) | `Figure_3.{pdf,png}` | `make_all_figures.py` ← `R3_pair_transformations/transformation_preference.ipynb` (`figure4_candidates.tsv`; per-cohort effects cached in `figure3_effects.csv`) | `manuscript/Figure4.png` |
| Figure 4 (non-CRC cohorts) | `Figure_4.{pdf,png}` | `make_all_figures.py` ← `R4_broader_application/analyse_synergy_gains.ipynb` (`synergy_by_cohort.tsv`) | `manuscript/Figure5.png` (without the two arrows of panel a) |
| Table 1 (17 CRC cohorts, sample counts) | manuscript text | `M1_data_preparation/prep_data.ipynb` (class counts per cohort) | — |
| Table 2 (non-CRC cohorts) | manuscript text | `R4_broader_application/analyse_synergy_gains.ipynb` | `manuscript/Table1_validation_cohorts.tsv` |
| Supp. Figure 1 (method coverage) | `Supplementary_Information.pdf` | `R5_reviewer_checks/run_method_coverage.py` | `manuscript/SuppFigure3_method_coverage.png` |
| Supp. Figure 2 (feature-only vs synthetic) | `Supplementary_Information.pdf` | `R2_classification_performance/summarize_rf_results.ipynb` | `manuscript/SuppFigure1_feature_only_vs_synthetic.png` |
| Supp. Figure 3 (AUC per test cohort) | `Supplementary_Information.pdf` | `R2_classification_performance/summarize_rf_results.ipynb` | `manuscript/SuppFigure2_auc_per_cohort.png` |
| Supp. Tables 1–9 | `Supplementary_Tables.xlsx` | `manuscript/build_supplementary_files.py` (sheet per table, from the TSVs listed in its `CAPTIONS`) | `manuscript/supplementary/Supplementary_tables.xlsx` (identical cell contents) |

`manuscript/Figure1.png`, `Figure1b.*` and `Figure3.*` are earlier versions that are not in the revised paper.
The submitted figures are rebuilt in the journal style with

```bash
cd manuscript/final
python make_all_figures.py --out . \
    --fig1-data ../../R1_synergistic_pairs/figure1b_all_pairs/directed \
    --fig2-data ../../R2_classification_performance \
    --fig3-candidates ../../R3_pair_transformations/figure4_candidates.tsv --fig3-cache figure3_effects.csv \
    --fig4-data ../../R4_broader_application/synergy_by_cohort.tsv
```

(`figure1b_all_pairs/directed/` holds the large per-fold tables; see *Data and large files* below.)

## Layout

| Folder | Manuscript section | Code | Produces |
|---|---|---|---|
| `run_mdfs.R` | Methods *MDFS feature selection and information gain* | one MDFS run: 1D IG + p-values, 2D significance, Bonferroni pair threshold, synergy filter | pair table, 1D table |
| `M1_data_preparation/` | Methods *Data preparation* | `prep_data.ipynb` | `prepared_data/X_*, y_*` (17 CRC cohorts × taxa / functions / combined; git-ignored) |
| `M2_loco_crc_pipeline/` | Methods *Leave-one-cohort-out (LOCO) validation* | `framework.ipynb` (functions), `main_code.ipynb` (orchestration), `run_loco.py` (headless) | `results/<dt>/`: `loco_results` (AUC), `loco_accuracy`, `base_features`, `synthetic_features` (RF selection frequencies), `mdfs_ig_info` (IGs of the top-2 % pairs of every fold × run), `mdfs_1d` (1D IG + p-values of every fold × run) |
| `R1_synergistic_pairs/` | Results *Synergy-oriented framework identifies disease-informative feature pairs* | `compare_mdfs_igs.ipynb`; `figure1b_all_pairs/` (all directed pairs, every fold × run); `seed_stability/` (20-seed stability grid) | **Figure 1b** (all directed pairs), **Supp. Table 2**, **Supp. Table 1**, Wilcoxon test, `ig_gain_summary_<dt>.csv` |
| `R2_classification_performance/` | Results *Synergy-derived features retain benchmark-level classification performance* | `summarize_rf_results.ipynb`; `stable_pairs.ipynb`; `effect_sizes.ipynb` | AUC summary (`loco_performance_summary.tsv`), **Supp. Figures 2–3**, stable-pair table (`final_synthetic_features_table.csv`), **Figure 2**, **Supp. Table 4**, meta-analysis tables |
| `R3_pair_transformations/` | Results *Pair transformations reveal distinct forms of disease-associated dependence* | `transformation_preference.ipynb` | **Figure 3**, **Supp. Table 6**, LR / GM counts |
| `R4_broader_application/` | Results *Broader application identifies conserved synergistic markers in other conditions*; Methods *Application to additional disease cohorts* | `run_mdfs_1d_2d.py` (10 MDFS runs per cohort), `load_cohort_data.py` (synergy table), `analyse_synergy_gains.ipynb` | **Figure 4**, **Table 2**, **Supp. Table 7**, `synergy_by_cohort.tsv`, emergent synergies, top markers per cohort |
| `R5_reviewer_checks/` | Reviewer comments 1.2, 1.5, 1.7, 1.9, 1.10, 2.1, 2.2, 2.5 | `run_negative_controls.py`, `run_discretisation_sensitivity.py`, `run_runtime.py`, `run_method_coverage.py`; `additivity_and_collinearity.ipynb`, `sensitivity_training_cohorts.ipynb`, `reviewer_checks_summary.ipynb` | label-permutation control, discretisation sensitivity, runtime/memory, cohort-adjusted likelihood-ratio tests of the stable pairs, collinearity, LOCO rerun without the single-class / training-only cohorts (`M2_loco_crc_pipeline/results_sensitivity/`), feature-level coverage of MDFS-2D versus MDFS-1D / t-test / U-test / LASSO / ReliefF (Supp. Figure 1), Supp. Tables 3, 5, 8, 9 |
| `manuscript/` | — | — | one file per manuscript figure / table (written by the notebooks above, never by hand; working names, see the table above), `REVISION_CHANGES.md` (submitted vs revised numbers); `final/`: the files as submitted and `make_all_figures.py` |
| `archive/` | — | analyses not used in the manuscript (see `archive/README.md`) | — |
| `METHODS.md` | — | step-by-step description of the whole pipeline with all parameters | — |

## Discretisation

MDFS bins every feature at a random cut point. `run_mdfs.R` uses **one seed and one cut range
(`GetRange(n, dimensions = 2)`) for the 1D, 2D and pair steps**, so a feature is binned identically in
every step and `total_IG = base_IG_1D + IG_2D_added` is the joint IG of the pair under one binning.
(The submitted version drew separate cuts for the 1D and 2D steps; see `manuscript/REVISION_CHANGES.md`.)

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
(≈ 3–5 h, the three data types run in parallel; add `python3 M2_loco_crc_pipeline/run_loco.py <dt> --exclude-training-only`
for the sensitivity run), `R1`–`R4` the result notebooks (minutes each, except `R2` effect sizes ≈ 10–60 min), `R5` the
reviewer checks (scripts first, then notebooks). Notebooks are executed in place with `jupyter nbconvert --execute`, so the
stored outputs are those of the last full run. The all-pairs input of Figure 1b
(`R1_synergistic_pairs/figure1b_all_pairs/run_all_directed.py`, several hours) is not part of
`run_all.sh`; its per-fold tables are git-ignored (below).

## Data and large files (not in git)

The input profiles are not redistributed in this repository. To rerun the analysis from the start,
place them as follows (paths relative to the repository root):

| Input | Expected location | Content |
|---|---|---|
| CRC discovery arm (17 cohorts; Table 1) | `../data/{taxa,pathways,metadata}/<cohort>.tsv`, i.e. a `data/` folder next to the repository, file stems `Public_study__<Study>` / `This_study_cohort<k>__<Cohort>` | MetaPhlAn species and HUMAnN pathway tables and sample metadata of the CRC meta-analysis benchmark (Piccinno et al. 2025) |
| Non-CRC cohorts (Table 2) | `R4_broader_application/data/` (`merged_abundances.txt`, `merged_pathways.txt`, `merged_metadata.tsv`, `metadata.tsv`) | curatedMetagenomicData species and pathway tables with metadata |

Everything else is produced by the pipeline. The large intermediate tables are git-ignored:
`M1_data_preparation/prepared_data/` (from `M1`), the per-run MDFS tables `mdfs_ig_info_*` / `mdfs_1d_*` of
`M2` and `R4` (15–95 MB each), `R1_synergistic_pairs/figure1b_all_pairs/directed/` and
`R1_synergistic_pairs/seed_stability/pairs/`. All summary tables and figures derived from them are in git.
An archive of the input profiles and of these intermediate tables from the run reported in the paper is
available from the authors on request.

**What runs from a fresh clone without the input data:** the supplementary tables
(`python manuscript/build_supplementary_files.py`; identical to the submitted file), Figures 2–4
(`manuscript/final/make_all_figures.py --only 2,3,4`, see above; Figure 3 uses the cached effect sizes), and
reading every notebook with its stored outputs. Figure 1b and all `M1`–`R5` computations need the inputs.

`manuscript/docx_tools/` contains the scripts used to edit the manuscript Word files during the revision.
They point to the author's local folders and are not needed to reproduce any result.

## Software

R 4.5 with `MDFS` 1.5.5 and `data.table` (`install.packages(c("MDFS", "data.table"))`); Python 3.9+ with the
packages in `requirements-manuscript.txt` (`pip install -r requirements-manuscript.txt`), which pins the
versions used for the manuscript (`METHODS.md` §7). The submitted figures were rendered with Liberation Sans;
with another font installed first (e.g. Arial) `make_all_figures.py` gives the same figures with slightly
different text rendering.
