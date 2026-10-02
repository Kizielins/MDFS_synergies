# Revision rerun: old vs new results

Rerun completed 2026-09-29. Original results are untouched; all new outputs are in `*_v2` folders.

## What changed in the method

| | Submitted | Revision |
|---|---|---|
| MDFS runs per LOCO fold | 2 (seed **not passed** to MDFS → both runs identical) | **10 independent runs**, seeds 42–51 |
| Pairs kept per run | top 1% by total IG | **top 2%** by total IG |
| Consensus (CRC) | in both runs | in **≥ 6 of 10** runs |
| Per-run diagnostic RF | trained each run (output unused) | removed |
| Final RF (AUC, selection frequencies) | unchanged | unchanged |
| Validation-arm runs | 3 (seeds 12–14), pair kept if found in **any** run | **10 runs** (seeds 12–21), pair kept if found in **≥ 6 of 10** |
| Feature / sample order | Python set (varied between sessions) | sorted (reproducible) |

Justification (Reviewer 1, comment 6): `../../seed_stability_analysis/consensus_grid.tsv`. With 10 runs and ≥ 6, the consensus sets built from two disjoint seed blocks agree with Jaccard 0.94 (taxa) and 0.96 (combined) at the 2% cutoff; with 2 runs they agreed at 0.53 and 0.56. Functions reach only ~0.35–0.41 at any cutoff (report as a limitation).

## CRC discovery arm

### Classification (Results §2, Abstract, Discussion; Supp. Fig. 1–2)
Mean AUC across the 14 test cohorts (weighted mean in brackets):

| | Raw features | Synthetic only | Raw + synthetic |
|---|---|---|---|
| Taxa | 0.849 → **0.847** (0.818) | 0.837 → **0.836** (0.809) | 0.837 → **0.838** (0.809) |
| Functions | 0.679 → **0.674** (0.664) | 0.673 → **0.667** (0.654) | 0.675 → **0.669** (0.657) |
| Combined | 0.821 → **0.825** (0.801) | 0.833 → **0.831** (0.803) | 0.831 → **0.824** (0.799) |

All changes are < 0.01. Text updates:
- Taxonomy "AUC of 0.84" and "synthetic-only 0.84 vs benchmark 0.85": unchanged.
- Combined: "synthetic and combined models on par (AUC = 0.83), marginally above raw (0.82)" → synthetic **0.83**, raw **0.83** (0.825), combined **0.82**. Suggest: "all three configurations reached AUC 0.82–0.83".
- Functions "AUC = 0.67–0.68" → **0.67** (all three).

Consensus pairs per fold: taxa 600 (533–1,104), functions 503 (326–692), combined 3,436 (3,087–4,494).
Methods text "~25,000–30,000 candidate pairs to ~250–300 per run" → per run the top 2% is ~550–1,100 pairs for taxa (from ~27,000–55,000 candidates).

### Figure 1b and IG comparison
Regenerated: `outputs/compare_mdfs_igs/figures/`. Wilcoxon 2D vs 1D IG: still p ≪ 0.001 for all three data types.

### Top synergistic pairs (Results §1, Discussion, Supp. Table 1)
New Supp. Table 1 source: `../results_<dt>_v2/ig_gain_summary_<dt>.csv` (`pair_ig_count` now out of **140** = 14 folds × 10 runs).

| Claim | Old | New |
|---|---|---|
| F. nucleatum in all top-ranking pairs | 10/10 top pairs by raw gain | 10/10 (unchanged) |
| F. nucleatum + P. stomatis | 26/28 | 135/140 |
| F. nucleatum + P. micra | 26/28 | 136/140 |
| F. nucleatum + D. pneumosintes | 26/28 | 133/140 |
| F. nucleatum + A. obesiensis | 26/28 | 130/140 |
| **P. micra + D. pneumosintes, "only pair at 28/28"** | **26/28 — the claim was incorrect in the submitted version**; the pair at 28/28 was D. pneumosintes + P. stomatis | **140/140** — at maximum, together with **P. micra + P. stomatis** (also 140/140) |

Suggested text: "P. micra + D. pneumosintes and P. micra + P. stomatis were the only pairs detected in all 140 fold–run combinations (14 folds × 10 MDFS runs)." (Results and Discussion; Reviewer 2 comment 6(a) response.)

### Stable pairs (≥ 80% of folds) — Results §2–3, Supp. Table 2
Source: `../downstream_analyses_v2/final_synthetic_features_table.csv`

| Claim | Old | New |
|---|---|---|
| Stable pairs | 62 | **211** |
| From combined model | 58 | **160** (taxa 34, functions 17) |
| ≥ 1 constituent at 0% individual frequency | 55 | **184** |
| P. micra raw frequency (combined) | 0% | **0%** (unchanged) |
| P. micra in stable combined pairs | 33 ("23 of 33 as base") | **19** (base/contributor split to recompute) |
| F. nucleatum raw frequency (combined) | 50% | **50%** (unchanged) |
| F. nucleatum stable combined partners | 16–17 | **114** |

More stable pairs is expected: the 2% cutoff gives the classifier ~2× more pairs per fold. The hub narrative shifts: **F. nucleatum becomes the dominant hub** (114 pairs); P. micra remains a 0%-individual hub captured through pairs (19).

### Effect-size meta-analysis — Results §3, Figure 3b
Source: `../downstream_analyses_v2/meta_results_all_pairs.csv`

| Claim | Old | New |
|---|---|---|
| Pairs meta-significant (Stouffer, BH q < 0.05) | 62/62 | **211/211** |
| Cohorts per pair | 15 | 15 |
| Consistent direction in ≥ 14/15 cohorts | 60 (27 at 15/15, 33 at 14/15; 2 at 13/15) | **201** (118 at 15/15, 83 at 14/15; 9 at 13/15, 1 at 12/15) |
| Cohort–pair tests significant (BH within pair) | 80.3% | **78.4%** (2,481/3,165) |

### Figure 3
`../downstream_analyses_v2/figure3_hub_representative_pairs.png`. Representative pairs now: **P. micra (base, 0%) + P. stomatis** (100% pair frequency, +100 pp) and **P. stomatis (base) + F. nucleatum (contributor, 50%)** (+75 pp); both significant in 14/15 cohorts. Legend text naming the representative pairs needs updating.

### Transformations — Results §4, Figure 4, Supp. Table 3
| Claim | Old | New |
|---|---|---|
| GM-based stable pairs | 41/62 (66%) | **127/211 (60%)** |
| Blautia wexlerae + P. micra, LR only, 100% | LR 100% | LR **36%**, **no longer stable** — replace example |
| P163-PWY + REDCITCYC, GM only | stable | still stable |
| P. micra + Eikenella corrodens, both forms | GM 100% | GM 86%, still stable |

Figure 4: `../downstream_analyses_v2/figure4_lr_vs_gm.png`. Supp. Table 3: `../downstream_analyses_v2/supplementary_table_3_transformation_preference.tsv` (746 pairs; was 566).

Not regenerated: `figure_synthetic_barplots.png` (one panel per stable pair; 211 panels exceed matplotlib's size limit). It is an overview, not a manuscript figure.

## Validation arm (non-CRC cohorts)
Source: `../apply_framework/synergy_by_cohort_v2.tsv`, figures in `../apply_framework/downstream_v2/`.

| Claim | Old | New |
|---|---|---|
| Unique synergistic pairs | 133,747 | **74,725** |
| ≥ 100% gain in ≥ 1 cohort | 112,888 (84.4%) | **69,573 (93.1%)** |
| Cohorts with significant pairs | 9 | **7** (NagySzakalD_2017, QinJ_2012 → 0; ZhuF_2020 → 1) |
| Pairs per cohort | QinN 72,385; HallAB 40,934; NagySzakalD 2; QinJ 8 | QinN 50,495; HallAB 14,916; JieZ 10,023; MetaCardis 3,729; RubelMA 447; NielsenHB 313; ZhuF 1 |
| QinN / JieZ with ≥ 100% gain | 99.8% / >99.9% | 99.8% / 100% |
| MetaCardis with ≥ 100% gain | 44.9% | 46.0% |
| Emergent synergies | 28,254 pairs; HallAB, MetaCardis, RubelMA | **5,171 pairs; HallAB, MetaCardis** |
| E. ramulus + R. torques (HallAB) | 0.58 → 55.6 (95-fold) | 1.31 → 52.9 (**~40-fold**), 9/10 runs |
| C. amalonaticus + S. anginosus (QinN) | 0.23 → 88.9 (~400-fold) | 0.23 → 90.2 (~390-fold), 9/10 runs |
| MetaCardis top pair | E. coli + L-Asp/L-Asn biosynthesis (61.4) | **E. coli + tRNA charging (61.5)** |
| **Gemmiger formicilis + adenine/adenosine salvage III** | 5 cohorts | **0 cohorts** (1–4 of 10 runs per cohort) |
| Conserved IBD pairs (HallAB & NielsenHB) | 176 (124 taxa–taxa, 52 taxa–function) | **11** |
| A. putredinis + G. pamelaeae ("most potent IBD") | present | **absent** (was in 1 run per cohort) |
| Max cohorts sharing a pair | 5 | **3** (94 pairs in 3 cohorts; 5,105 in ≥ 2) |

Most-recurrent pairs now (3 cohorts: HallAB, JieZ, QinN) pair an oral Streptococcus (S. vestibularis, S. mitis, S. anginosus group, S. oralis) with a gut commensal.

Sections to rewrite: Results "Broader application…" (counts, emergent, examples, Gemmiger, IBD), Discussion Gemmiger paragraph, Abstract ("over 130,000 synergistic pairs"), Supp. Table 4, Figure 5 (consider dropping cohorts with < 10 pairs), Reviewer 1 comment 12 response.

## Methods / response-letter text to change
- LOCO paragraph: "MDFS was run 2 times … (42, 43) … empirically the two runs were highly stable" → 10 independent runs (seeds 42–51), top 2% per run, consensus ≥ 6/10; cite the stability analysis.
- Results §1: "run twice per fold … top 1% of IG in both runs" → same update.
- Delete "RF consensus across runs was defined analogously to MDFS…".
- "random_state matching the run seed (42 or 43)" → random_state = 42.
- Validation arm: "three random seeds (12, 13, 14); pair-level IG values were averaged across seeds" → 10 seeds (12–21); pairs retained if significant in ≥ 6/10 runs; IGs averaged over those runs.
- Figure 1b legend: "top 1% of the features" → top 2%.
- Reviewer 1 comment 6: new answer + stability table as a supplementary table.
- Reviewer 2 comment 6(a): "28/28 … 14 folds × 2 MDFS runs" → 140/140, 14 × 10.
- Reviewer 2 minor 4: the tool is on `main`; the manuscript analysis code is on the `manuscript-analysis` branch.

## Where everything is
| Output | Location (relative to `framework_repository/`) |
|---|---|
| LOCO results | `results_taxa_v2/`, `results_functions_v2/`, `results_combined_v2/` |
| Fig. 1b, Supp. Table 1 | `rerun_v2/outputs/compare_mdfs_igs/figures/`, `results_*_v2/ig_gain_summary_*.csv` |
| Supp. Fig. 1–2 | `rerun_v2/outputs/summarize_rf_results/figures/` |
| Fig. 3, Fig. 4, Supp. Tables 2–3, meta-analysis | `downstream_analyses_v2/` |
| Incomplete pairings | `rerun_v2/outputs/analyse_synthetic_combinations/` |
| Fig. 5, Supp. Table 4, emergent synergies | `apply_framework/downstream_v2/`, `apply_framework/synergy_by_cohort_v2.tsv` |
| Run logs | `rerun_v2/*.log` |
| Seed-stability analysis | `../seed_stability_analysis/` |
