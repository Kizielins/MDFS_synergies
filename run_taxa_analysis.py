#!/usr/bin/env python3
"""
Script to run LOCO analysis for taxa data with feature intersection approach.
This extracts and runs the main function from main_code.ipynb
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split
import warnings
import subprocess
import os
import matplotlib.pyplot as plt
from collections import Counter

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# --- CONFIGURATION ---
ALL_COHORT_STEMS = [
    'Public_study__FengQ_2015', 'Public_study__GuptaA_2019', 'Public_study__LiuNN_2022',
    'Public_study__ThomasAM_2018b', 'Public_study__VogtmannE_2016', 'Public_study__WirbelJ_2018',
    'Public_study__YachidaS_2019', 'Public_study__YangJ_2020', 'Public_study__YangY_2021',
    'Public_study__YuJ_2015', 'Public_study__ZellerG_2014', 'This_study_cohort1__ONCOBIOME_AtezoTRIBE',
    'This_study_cohort2__ONCOBIOME_COLOBIOME', 'This_study_cohort3__ONCOBIOME_IIGM_CZ',
    'This_study_cohort4__ONCOBIOME_IIGM_IT', 'This_study_cohort5__NSHII', 'This_study_cohort6__IIGM_TU'
]
TRAINING_ONLY_COHORT = 'This_study_cohort5__NSHII'

# --- MODEL AND FEATURE SELECTION PARAMETERS ---
TOP_K_VALUES = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
MODEL_PARAMS = {
    'n_estimators': 1000, 'max_features': 'sqrt', 'min_samples_leaf': 1,
    'class_weight': 'balanced', 'random_state': 42, 'n_jobs': -1
}

SYNTHETIC_PREFIXES = ('LR_', 'GM_')

# Import functions from the notebook (we'll need to define them here or import)
# For now, let's just run the main analysis
if __name__ == '__main__':
    print("Running taxa analysis with feature intersection...")
    print("Note: This requires the functions from main_code.ipynb to be available.")
    print("Please run the notebook cells in order, or import the functions first.")
