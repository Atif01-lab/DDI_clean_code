# outputs/

Created and populated when the notebooks are run (location set by `OUTPUT_DIR` in `analysis/paths.py`). Nothing in this folder is committed to the repository.

| Folder | Contents | Shareable? |
|---|---|---|
| `derived/` | Cohort, event, covariate and polypharmacy files | **No - patient-level data** |
| `results/` | Pooled hazard-ratio summaries and reference strata (main analyses) | Aggregate only, subject to CPRD small-number rules |
| `supplementary/` | Bleeding-site and DOAC-type results | Aggregate only, subject to CPRD small-number rules |
| `figures/main/`, `figures/supplementary/` | Forest plots and other figures | Yes, once checked |
| `sensitivity/` | Results of the four sensitivity analyses | Aggregate only, subject to CPRD small-number rules |
