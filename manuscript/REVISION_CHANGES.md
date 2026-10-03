# Revision: submitted version vs this revision

Generated from the final output files (`manuscript/collect_numbers.py` values). The submitted numbers are those of the manuscript sent to npj Biofilms and Microbiomes; 'September rerun' is the intermediate 10-run rerun whose numbers were used in `MS_final_reviewed.docx`.

## What changed in the method

| | Submitted | This revision |
|---|---|---|
| MDFS runs per LOCO fold | 2, seed not passed to MDFS (runs identical) | 10 independent runs, seeds 42-51 |
| Pairs kept per run / consensus | top 1 %, in both runs | top 2 %, in >= 6 of 10 runs (seed-stability analysis, Supp. Table 5) |
| Discretisation | 1D and 2D steps used different random cut points | one seed and one cut range (GetRange(n, 2)) for the 1D, 2D and pair steps |
| Per-run diagnostic RF | trained, output unused | removed; one RF per fold and configuration (random_state 42) |
| Validation arm | 3 seeds, pair kept if found in any run; ordered (directed) pair counts | 10 seeds (12-21), kept if in >= 6 of 10 runs; unordered pairs, orientation recorded |
| 1D p-values | not saved | saved for every fold and run (`mdfs_1d_<dt>.csv`) |
| New analyses | - | negative controls, discretisation sensitivity, runtime, interaction information, cohort-adjusted LRT, collinearity, training-composition sensitivity (R5_reviewer_checks/) |

## CRC discovery arm

| Claim | Submitted | This revision |
|---|---|---|
| Mean AUC taxa (raw / synthetic / both) | 0.849 / 0.837 / 0.837 | 0.847 / 0.839 / 0.837 |
| Mean AUC functions | 0.679 / 0.673 / 0.675 | 0.674 / 0.667 / 0.672 |
| Mean AUC combined | 0.821 / 0.833 / 0.831 | 0.825 / 0.824 / 0.827 |
| Stable pairs (combined / taxa / functions) | 62 (58 / - / -) | 235 (169 / 40 / 26) |
| >= 1 constituent at 0 % individual frequency | 55 | 207 |
| P. micra stable combined pairs (as base) | 33 (23) | 18 (13) |
| F. nucleatum stable combined pairs (as contributor) | 16-17 | 125 (123) |
| Pairs meta-significant / consistent direction >= 14/15 | 62/62 / 60 | 235/235 / 225 |
| Cohort-pair tests significant | 80.3 % | 78.0 % |
| GM-based stable pairs (GM only + both) | 41/62 (66 %) | 173/235 (74 %) |
| Most consistent pairs | P. micra + D. pneumosintes 'only pair at 28/28' (incorrect) | taxa: F. nucleatum + P. micra (140/140), P. micra + P. stomatis (140/140); combined: G. adiacens + P. stomatis (140/140), A. obesiensis + G. morbillorum (140/140) |
| Seed stability (Jaccard, 10 runs >= 6 vs 2 runs) | - | taxa 0.95 vs 0.47; combined 0.98 vs 0.83; functions 0.33 vs 0.13 |
| P. micra 1D-significant (combined / taxa folds) | 'all folds' | 14/14 / 14/14 |
| Interaction information > 0 (consensus pairs) | - | taxa 49.7 % of 77,060 pair-fold entries (median II -0.0 bits, -0 % of the pair total IG); combined 39.2 % of 506,761 pair-fold entries (median II -1.4 bits, -0 % of the pair total IG) |
| Cohort-adjusted LRT (stable pairs beyond cohort / interaction term) | - | 235/235 / 123/235 |
| Sensitivity run (no single-class / training-only cohorts): taxa AUC | - | 0.845 / 0.836 / 0.831 (main 0.847 / 0.839 / 0.837); stable pairs shared: 169 / 135 / 61 (0.25) |
| Negative controls (max pairs: label permutation / within-class shuffle) | - | 31,751 vs 0.0 (max 0); 31,751 vs 30645.8 (max 30671) |

## Validation arm (non-CRC cohorts)

| Claim | Submitted | This revision |
|---|---|---|
| Cohorts analysed / with pairs | '13' / 9 | 13 / 8 |
| Unique pairs | 133,747 | 80,790 |
| >= 100 % gain in >= 1 cohort | 112,888 (84.4 %) | 76,828 (95.1 %) |
| Pairs per cohort | QinN 72,385; HallAB 40,934 | QinN_2014 55,853; HallAB_2017 12,077; JieZ_2017 11,137; MetaCardis_2020_a 7,585; NielsenHB_2014 313; RubelMA_2020 156; NagySzakalD_2017 2; ZhuF_2020 1 |
| Emergent pairs | 28,254 in 3 cohorts | 3,808 in 2 (HallAB_2017, MetaCardis_2020_a) |
| E. ramulus + R. torques (HallAB) | 0.58 -> 55.6 (95-fold) | 0.58 -> 54.9 (94-fold, 8/10 runs) |
| C. amalonaticus + S. anginosus (QinN) | 0.23 -> 88.6 | 0.23 -> 88.2 (10/10 runs) |
| MetaCardis top pair | E. coli + L-Asp/L-Asn biosynthesis (61.4) | ASPASN-PWY:_superpathway_of_L-aspartate_and_L-asparagine_biosynthesis + PWY-6609:_adenine_and_adenosine_salvage_III (61.0) |
| Gemmiger formicilis + PWY-6609 | 5 cohorts | 1 cohort(s) |
| Max cohorts sharing a pair | 5 | 3 (102 pairs; 89 % with Streptococcus) |
| Conserved IBD pairs | 176 (124 + 52) | 14 (10 + 4); same orientation in both cohorts: 0 |
| Table 1 | disease / control columns swapped | corrected (manuscript/Table1_validation_cohorts.tsv) |

## Where everything is

See `README.md` (section -> folder -> file) and `manuscript/RESULTS_NUMBERS.md` (every quoted number with its source). The AI-edited manuscript and reviewer-response drafts are outside the repository (NPJ_review/*AI_revision*.docx).
