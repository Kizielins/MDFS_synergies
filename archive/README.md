# Archive

Analyses that were explored during the project but are **not used by any figure, table or number in the manuscript**. Kept for reference only; paths were updated to the final repository layout and stored outputs were cleared.

| Item | What it is |
|---|---|
| `analyse_synthetic_combinations.ipynb` | "Incomplete pairings": for every stable pair, which of the two encodings (LR, GM) the RF selected, with forest plots of per-cohort effect sizes of the chosen and the never-chosen encoding. The LR/GM comparison used in the manuscript is `R3_pair_transformations/transformation_preference.ipynb`. |
| `investigate_synergies/` | Hand-curated check (`pairs.csv`) of a few pairs: 1D and 2D IGs of each listed pair on the pooled CRC data (`investigate_synergies.r`, outputs `ig_synergy_<data_type>.csv`). No random seed is set, so the IGs are not reproducible bit for bit. Run from this folder with `Rscript investigate_synergies.r --prepared_data_dir ../../M1_data_preparation/prepared_data --data_type combined --pairs pairs.csv --out ig_synergy_combined.csv`. |
