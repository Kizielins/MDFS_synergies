#!/usr/bin/env python3
"""Assemble the supplementary files from the manuscript/ outputs.

  Supplementary_tables.xlsx   one sheet per Supplementary Table (1-5), caption in row 1, data from the TSVs
  Supplementary_Figures.docx  copy of NPJ_review/Supplementary_Figures.docx with the two images replaced by
                              SuppFigure1_*.png and SuppFigure2_*.png (extent rescaled to the new aspect ratio)
Nothing is written into NPJ_review/; outputs go to manuscript/supplementary/.
"""
import os
import re
import shutil
import struct
import zipfile
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "supplementary")
os.makedirs(OUT, exist_ok=True)
NPJ = "/Users/kizie/Library/Mobile Documents/com~apple~CloudDocs/OneDrive-UJ/Sonata2/NPJ_review"

CAPTIONS = {
 1: ("SuppTable2_top30_synergistic_pairs.tsv", "Supplementary Table 2. The 30 directed synergistic feature pairs with the highest mean total 2D information gain (IG) across the taxa-only, functions-only and combined analyses. base_IG_1D_mean: 1D IG of the base feature; total_IG_2D_mean: joint IG of the pair (base 1D IG + IG added by the contributing feature); ig_difference_raw / ig_difference_pct: absolute and relative gain over the base feature; pair_ig_count: number of fold-run combinations (of 140 = 14 LOCO folds x 10 MDFS runs) in which the pair was among the top 2 % of pairs."),
 2: ("SuppTable4_stable_pairs.tsv", "Supplementary Table 4. Stable synthetic feature pairs in the CRC data: pairs whose synthetic feature was selected by the final random forest in >= 80 % of the 14 LOCO folds and more often than either constituent (Selection % Gain = pair frequency minus the mean frequency of its constituents). Post-selection effect-size summary per pair: transformation evaluated in most cohorts, number of two-class cohorts, Benjamini-Hochberg q of the Stouffer meta p-value across cohorts, median rank-biserial correlation (positive = higher in CRC) and the number of cohorts with positive / negative direction."),
 3: ("SuppTable6_transformation_preference.tsv", "Supplementary Table 6. Transformation preference of synthetic feature pairs: for every pair selected in >= 50 % of LOCO folds as a log-ratio (LR) and/or a geometric mean (GM), the encodings that reached this threshold (Types_Present) and the maximum selection frequency."),
 4: ("SuppTable7_conserved_pairs.tsv", "Supplementary Table 7. Candidate conserved markers: synergistic pairs identified independently in >= 2 cohorts with the same diagnosis (only inflammatory bowel disease, HallAB_2017 and NielsenHB_2014, met this criterion). Pairs were matched irrespective of orientation; orientation_per_cohort gives the (base -> contributing) direction reported in each cohort and same_orientation_in_all_cohorts whether it agrees. IG values are means over the cohorts."),
 6: ("../R5_reviewer_checks/negative_controls.tsv", "Supplementary Table 8. Negative controls for the pair search (LOCO fold-0 training sets, seed 42): number of 1D- and 2D-significant features and of synergistic pairs returned by run_mdfs.R with the real labels, with labels permuted across samples (five permutations) and with every feature shuffled independently within each class (five permutations; preserves every univariate signal and destroys all dependence between features)."),
 7: ("../R5_reviewer_checks/discretisation_sensitivity.tsv", "Supplementary Table 9. Sensitivity of the pair search to the MDFS discretisation (LOCO fold-0 training set, taxonomic features): number of synergistic pairs and overlap with the default setting (Jaccard index of all pairs and of the top-2 % pairs; Spearman correlation of the joint IG of shared pairs) when the admissible range of cut points is widened to the one-dimensional default (range_1D) or halved (range_half), when three bins per feature are used (divisions_2) and, for reference, when only the random seed is changed (default_seed43)."),
 8: ("../R5_reviewer_checks/sensitivity_auc.tsv", "Supplementary Table 3. Sensitivity of the LOCO results to the training-set composition: mean and sample-size-weighted AUC of the three feature configurations with the main training sets and with the three cohorts that are never held out (AtezoTRIBE, COLOBIOME, NSHII) removed from every training set. The stable pairs shared by both runs are listed in R5_reviewer_checks/sensitivity_shared_stable_pairs.tsv."),
 9: ("../R5_reviewer_checks/stable_pairs_lrt.tsv", "Supplementary Table 5. Cohort-adjusted logistic regression for the stable pairs (pooled samples of the 15 two-class CRC cohorts, log10 relative abundances, cohort as fixed effect): likelihood-ratio p-values of the additive model against the cohort-only model (p_pair_beyond_cohort) and of the model with a multiplicative interaction term against the additive model (p_interaction), with Benjamini-Hochberg q-values."),
 5: ("SuppTable1_seed_stability.tsv", "Supplementary Table 1. Seed stability of the MDFS pair selection. For each data type, number of MDFS runs N and consensus threshold k (pair kept if in the per-run top 2 % of >= k of N runs): mean number of consensus pairs per LOCO fold and the Jaccard agreement (mean and minimum over the 14 folds) between consensus sets built from two disjoint blocks of seeds (42-51 and 52-61). The pipeline uses N = 10, k = 6."),
}
wb = Workbook(); wb.remove(wb.active)
for f, cap in sorted(CAPTIONS.values(), key=lambda fc: int(re.search(r"Supplementary Table (\d+)", fc[1]).group(1))):
    n = int(re.search(r"Supplementary Table (\d+)", cap).group(1))
    p = os.path.normpath(os.path.join(HERE, f))
    if not os.path.exists(p):
        print("missing", f); continue
    df = pd.read_csv(p, sep="\t")
    ws = wb.create_sheet(f"Supplementary Table {n}")
    ws.cell(1, 1, cap).font = Font(bold=True); ws.cell(1, 1).alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=max(2, min(len(df.columns), 8)))
    ws.row_dimensions[1].height = 90
    for j, c in enumerate(df.columns, 1):
        ws.cell(3, j, c).font = Font(bold=True)
    for i, row in enumerate(df.itertuples(index=False), 4):
        for j, v in enumerate(row, 1):
            ws.cell(i, j, None if pd.isna(v) else (float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else str(v)))
    for j, c in enumerate(df.columns, 1):
        ws.column_dimensions[ws.cell(3, j).column_letter].width = min(60, max(12, len(str(c)) + 2))
wb.save(os.path.join(OUT, "Supplementary_tables.xlsx")); print("wrote Supplementary_tables.xlsx")


def png_size(path):
    with open(path, "rb") as fh:
        head = fh.read(24)
    return struct.unpack(">II", head[16:24])


src = os.path.join(NPJ, "Supplementary_Figures.docx")
dst = os.path.join(OUT, "Supplementary_Figures.docx")
imgs = {"word/media/image1.png": os.path.join(HERE, "SuppFigure1_feature_only_vs_synthetic.png"),
        "word/media/image2.png": os.path.join(HERE, "SuppFigure2_auc_per_cohort.png")}
zin = zipfile.ZipFile(src)
doc = zin.read("word/document.xml").decode("utf8")
rels = zin.read("word/_rels/document.xml.rels").decode("utf8")
rid = {m.group(2): m.group(1) for m in re.finditer(r'Id="(rId\d+)"[^>]*Target="(media/image\d+\.png)"', rels)}
for media, new in imgs.items():
    r = rid[media.split("word/")[1]]
    w, h = png_size(new)
    # rescale the drawing extent that embeds this image: keep width, set height from the new aspect ratio
    def fix(m):
        block = m.group(0)
        ext = re.search(r'<wp:extent cx="(\d+)" cy="(\d+)"/>', block)
        cx = int(ext.group(1)); cy = int(round(cx * h / w))
        block = block.replace(ext.group(0), f'<wp:extent cx="{cx}" cy="{cy}"/>')
        block = re.sub(r'<a:ext cx="(\d+)" cy="\d+"/>', lambda mm: f'<a:ext cx="{mm.group(1)}" cy="{int(round(int(mm.group(1)) * h / w))}"/>', block)
        return block
    doc = re.sub(r'<w:drawing>(?:(?!</w:drawing>).)*r:embed="%s"(?:(?!</w:drawing>).)*</w:drawing>' % r, fix, doc, flags=re.S)
with zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        if item.filename == "word/document.xml":
            zout.writestr(item, doc.encode("utf8"))
        elif item.filename in imgs:
            zout.write(imgs[item.filename], item.filename)
        else:
            zout.writestr(item, zin.read(item.filename))
print("wrote Supplementary_Figures.docx (images replaced)")
