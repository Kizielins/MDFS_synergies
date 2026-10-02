#!/usr/bin/env python3
"""
Headless run of the CRC LOCO pipeline for one data type (Methods: LOCO validation).

Executes the code cells of framework.ipynb and main_code.ipynb (except the
notebook-loading and run-all cells) and calls main(); results are written to
results/<data_type>/. Each data type runs in its own working directory
(work_<data_type>/, git-ignored) so the MDFS temporary files of parallel runs
cannot collide.

Usage: python run_loco.py <taxa|functions|combined> [--exclude-training-only]
  --exclude-training-only  sensitivity run (Reviewer 1, comment 5): the three cohorts that are never held
                           out (AtezoTRIBE, COLOBIOME: CRC only; NSHII: training-only) are also dropped from
                           every training set; results go to results_sensitivity/<data_type>/.
"""
import json
import os
import shutil
import sys

HERE = os.path.abspath(os.path.dirname(__file__))
REPO = os.path.dirname(HERE)
data_type = sys.argv[1]
assert data_type in ("taxa", "functions", "combined")
SENS = "--exclude-training-only" in sys.argv
RESULTS_DIR = "results_sensitivity" if SENS else "results"


def code_cells(path):
    nb = json.load(open(path, encoding="utf-8"))
    return ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]


workdir = os.path.join(HERE, f"work_{data_type}" + ("_sensitivity" if SENS else ""))
os.makedirs(workdir, exist_ok=True)
shutil.copy(os.path.join(REPO, "run_mdfs.R"), workdir)
os.chdir(workdir)

g = {"__name__": "run_loco"}
main_cells = code_cells(os.path.join(HERE, "main_code.ipynb"))
exec(main_cells[0], g)                                   # imports + cohort config
for src in code_cells(os.path.join(HERE, "framework.ipynb")):
    exec(src, g)                                         # framework functions
exec(main_cells[2], g)                                   # load_prepared_data
exec(main_cells[3], g)                                   # main()
if SENS:
    g["EXCLUDE_FROM_TRAINING"] = ("This_study_cohort1__ONCOBIOME_AtezoTRIBE", "This_study_cohort2__ONCOBIOME_COLOBIOME", "This_study_cohort5__NSHII")
    print("SENSITIVITY RUN: excluded from every training set:", g["EXCLUDE_FROM_TRAINING"])

g["main"](
    prepared_data_dir=os.path.join(REPO, "M1_data_preparation", "prepared_data"),
    data_type=data_type,
    use_synthetic=True,
    scenario="All Features + Synthetic",
    output_dir=os.path.join(HERE, RESULTS_DIR, data_type),
)
print("RUN COMPLETE", data_type)
