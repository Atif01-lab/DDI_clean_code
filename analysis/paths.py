"""
Central file locations for the DOAC drug-drug interaction (DDI) analysis.

Every notebook in this folder starts with ``from paths import *``. Edit the two
root folders below, or set the environment variables ``DDI_DATA_DIR`` and
``DDI_OUTPUT_DIR`` before starting Jupyter, so that no notebook needs editing.

DATA_DIR    - restricted CPRD Aurum extract, lookups and study code lists
              (NOT distributed with this repository; see data/README.md)
OUTPUT_DIR  - everything the notebooks write. Files in OUTPUT_DIR/derived are
              patient-level and must never be shared or committed.
"""
import os

REPO_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA_DIR = os.environ.get("DDI_DATA_DIR", os.path.join(REPO_DIR, "data"))
OUTPUT_DIR = os.environ.get("DDI_OUTPUT_DIR", os.path.join(REPO_DIR, "outputs"))

# ---- Inputs: CPRD Aurum extract and lookups (restricted) ---------------------
DRUG_ISSUE_ZIP = os.path.join(DATA_DIR, "DrugIssue.zip")      # DrugIssue files (tab-delimited .txt)
OBSERVATION_ZIP = os.path.join(DATA_DIR, "obs_data.zip")      # Observation files (tab-delimited .txt)
PATIENT_FILE = os.path.join(DATA_DIR, "pats_data", "pats_data.txt")  # Patient file (comma-delimited)
PRODUCT_DICT = os.path.join(DATA_DIR, "prod_dict.txt")        # CPRD Aurum product dictionary
LOOKUP_DIR = os.path.join(DATA_DIR, "Lookup")                 # CPRD lookups / medical dictionary
DUPLICATE_PRACTICES_FILE = os.path.join(LOOKUP_DIR, "duplicate_practice_ids.txt")

# ---- Inputs: study code lists (see supplementary tables S1, S3-S5) -----------
CODE_LIST_DIR = os.path.join(DATA_DIR, "code_list")

# ---- Outputs -----------------------------------------------------------------
DERIVED_DIR = os.path.join(OUTPUT_DIR, "derived")             # patient-level - never share
RESULTS_DIR = os.path.join(OUTPUT_DIR, "results")             # aggregate model results
SUPP_DIR = os.path.join(OUTPUT_DIR, "supplementary")          # aggregate supplementary results
FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures", "main")
SUPP_FIGURES_DIR = os.path.join(OUTPUT_DIR, "figures", "supplementary")
SENS_WARFARIN_DIR = os.path.join(OUTPUT_DIR, "sensitivity", "warfarin_30d")
SENS_GRACE_DIR = os.path.join(OUTPUT_DIR, "sensitivity", "grace_period_30d")
SENS_COMPETING_DIR = os.path.join(OUTPUT_DIR, "sensitivity", "competing_risk")
SENS_LAG_DIR = os.path.join(OUTPUT_DIR, "sensitivity", "lag_exposure")

for _d in (DERIVED_DIR, RESULTS_DIR, SUPP_DIR, FIGURES_DIR, SUPP_FIGURES_DIR,
           SENS_WARFARIN_DIR, SENS_GRACE_DIR, SENS_COMPETING_DIR, SENS_LAG_DIR):
    os.makedirs(_d, exist_ok=True)

__all__ = [n for n in dir() if n.isupper()]
