# data/

**This folder is intentionally empty.** The study used Clinical Practice Research Datalink (CPRD) Aurum data, which are available only to approved researchers under a CPRD licence and an approved study protocol; they cannot be shared. Researchers with access can request equivalent data from CPRD (https://www.cprd.com).

Point `DATA_DIR` in `analysis/paths.py` (or the `DDI_DATA_DIR` environment variable) at a folder with the layout below. File names are those expected by the code.

```
DATA_DIR/
├── DrugIssue.zip            CPRD Aurum DrugIssue files (tab-delimited .txt inside the archive)
├── obs_data.zip             CPRD Aurum Observation files (tab-delimited .txt inside the archive)
├── prod_dict.txt            CPRD Aurum product dictionary (tab-delimited)
├── pats_data/
│   └── pats_data.txt        CPRD Aurum Patient file (comma-delimited)
├── Lookup/
│   ├── med_dict.txt                 CPRD Aurum medical dictionary
│   ├── common_dosages.txt           CPRD common dosages lookup (dosageid -> daily_dose)
│   ├── EmisCodeCat.txt, ObsType.txt, NumUnit.txt, Gender.txt   CPRD lookups
│   ├── snomed_icd10_mappings.txt    SNOMED CT -> ICD-10 mapping (columns SNOMED, ICD10_CODE)
│   └── duplicate_practice_ids.txt   CPRD practice IDs flagged as duplicates (one per line)
└── code_list/               Study code lists (content reported in supplementary tables S1, S3-S5)
    ├── bleeding_codes.xlsx              MedCodeId, Bleeding Category
    ├── hospitalisation_codes.xlsx       MedCodeId
    ├── comorbidities_codes.xlsx         one sheet per covariate, column MedCodeId
    ├── DOAC_DrugInteractions_Table.xlsx Drug, Category, Action, PK/PD  (table S1)
    ├── DOAC_interacting_drugs.xlsx      ProdCodeId, DrugSubstanceName, matched_drugs, Term from EMIS  (table S3)
    └── CodeList/
        └── covariate_ethniciity_aurum.txt  medcodeid, eth5
```

Notes

- Fields used from the Patient file: `patid`, `pracid`, `gender`, `yob`, `regstartdate`, `regenddate`, `cprd_ddate`, `acceptable`.
- Fields used from DrugIssue: `patid`, `issuedate`, `prodcodeid`, `dosageid`, `quantity`, `duration`, `probobsid`.
- Fields used from Observation: `patid`, `obsdate`, `medcodeid`, `value`, `numunitid`, `obstypeid`.
- `prod_dict.txt` must include the CPRD columns `ProdCodeId`, `DrugSubstanceName` and `Term from EMIS`, plus a `drug_class` column used to identify proton pump inhibitors (value containing "Proton Pump Inhibitor").
- `comorbidities_codes.xlsx` sheets: hypertension, stroke_tia_se, diabetes, cancer, smoking, chd, alcohol, anemia, heart_failure, angina, myocardial_infarction, atrial_fibrillation, vte, liver, ckd, bmi, hemoglobin, platelet_count, cholestrol_total, cholestrol_ldl, cholestrol_hdl, systolic_bp, diastolic_bp, pulse_rate.
