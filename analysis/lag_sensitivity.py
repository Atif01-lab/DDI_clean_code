"""
Exposure-lag sensitivity analyses for the DOAC drug-drug interaction study.
==========================================================================

Purpose
-------
Assesses protopathic bias: if a patient is classified as exposed from the day an
interacting drug is issued, a bleed that was already developing may appear to be
associated with the drug. Exposure onset is therefore lagged, applying the same
lag length to every drug / drug class (Table S17 reports the 3-day lag).

What this module does
---------------------
It re-runs the whole exposure-classification pipeline with the interacting-drug
exposure lagged by L days, for L in {0, 3, 7, 14}, applying the SAME lag to
every drug / drug class. L = 0 reproduces the primary analysis and exists as a
self-check: if the lag-0 table does not match the published table, the rebuild
path is wrong and none of the lagged numbers should be believed.

Two lag modes are supported (`mode=` argument of `lag_episodes`):

  "delay_onset"  (default) - an episode [s, e) becomes [s + L, e). Exposure
        starts L days after the prescription; the end (including the grace
        period already applied by merge_prescriptions) is unchanged. Episodes
        shorter than L disappear entirely and that person-time reverts to the
        drug-free reference. This is the literal reading of "the drug takes
        time to work", and it is the version to report as the primary
        sensitivity analysis.

  "shift_window" - an episode [s, e) becomes [s + L, e + L), i.e. E(t - L).
        Total exposed person-time is preserved, so short courses are not
        destroyed, but the lag doubles as extra carry-over after the last dose.
        Useful as a second row in the supplement, especially for clarithromycin
        / fluconazole where "delay_onset" removes most of the exposed time.

IMPORTANT - what the lag does NOT change
----------------------------------------
* The cohort, follow-up time, events and censoring are untouched. Only the
  classification of person-time as exposed / unexposed changes.
* Because follow-up and events are unchanged, the multiple imputations built
  for the primary analysis are reused by default (`imputations=` argument).
  The imputation model does contain ever-exposed flags, which shift slightly
  under a lag, but the imputed variables (eGFR, BMI) are baseline covariates and
  the effect is negligible. Set `refit_mice=True` (and supply the MICE helpers)
  if imputation is to be repeated inside each lag.

Reporting
---------
`exposure_retention()` produces the table that MUST accompany any lagged
estimate: exposed events and exposed person-years retained at each lag. A lagged
hazard ratio computed on 4 remaining events is not evidence of "no effect"; it
is an uninformative analysis and has to be labelled as such.

Author note: every function here is pure pandas/numpy. The notebook's own
`build_exposure_timeline`, Cox fitters and summary builders are passed in, so
this module never duplicates the modelling code.
"""

from __future__ import annotations

import gc
import time

import numpy as np
import pandas as pd

__all__ = [
    "LAG_DAYS",
    "lag_episodes",
    "build_lagged_frame",
    "events_and_py",
    "exposure_retention",
    "run_lag_sensitivity",
    "combine_lag_summaries",
    "lag_comparison_wide",
]

# Lags applied universally to every drug / class. 0 = primary analysis (control).
LAG_DAYS = [0, 3, 7, 14]


# ---------------------------------------------------------------------------
# 1. Lagging the merged exposure episodes
# ---------------------------------------------------------------------------
def lag_episodes(presc_merged, lag_days, term_col="drug_name", mode="delay_onset",
                 label=""):
    """Apply a universal lag to already-merged exposure episodes.

    The lag is applied AFTER `merge_prescriptions`, so the 7-day (antibiotic /
    antifungal) and 14-day (all other classes) gaps that define a continuous
    exposure episode are unaffected - the lag shifts the episode that those
    rules produced, it does not change the rules.

    Parameters
    ----------
    presc_merged : DataFrame
        Output of `merge_prescriptions` with `study_start_date` /
        `study_end_date` attached: one row per continuous exposure episode with
        columns patid, `term_col`, issuedate, enddate, study_start_date,
        study_end_date.
    lag_days : int
        Universal lag L in days. 0 returns the episodes unchanged (bar the
        clipping, which is a no-op on already-clipped input).
    term_col : {"drug_name", "drug_class"}
        Which exposure the notebook models.
    mode : {"delay_onset", "shift_window"}
        See the module docstring.

    Returns
    -------
    DataFrame with the same columns, episodes with zero or negative length
    dropped, all dates clipped to the patient's DOAC follow-up window.
    """
    if mode not in ("delay_onset", "shift_window"):
        raise ValueError("mode must be 'delay_onset' or 'shift_window'")
    if lag_days < 0:
        raise ValueError("lag_days must be >= 0")

    df = presc_merged.copy()
    for c in ("issuedate", "enddate", "study_start_date", "study_end_date"):
        df[c] = pd.to_datetime(df[c])

    n_in = len(df)
    lag = pd.to_timedelta(int(lag_days), unit="D")

    df["issuedate"] = df["issuedate"] + lag
    if mode == "shift_window":
        df["enddate"] = df["enddate"] + lag

    # Re-clip to the DOAC window. Shifting can push an episode past study end.
    df["issuedate"] = df["issuedate"].clip(lower=df["study_start_date"],
                                           upper=df["study_end_date"])
    df["enddate"] = df["enddate"].clip(lower=df["study_start_date"],
                                       upper=df["study_end_date"])

    keep = df["issuedate"] < df["enddate"]
    dropped = int((~keep).sum())
    df = df.loc[keep].reset_index(drop=True)

    exposed_days = (df["enddate"] - df["issuedate"]).dt.total_seconds().sum() / 86400.0
    print(f"[{label}] lag {lag_days:>2}d ({mode}): {len(df):,} episodes kept of "
          f"{n_in:,} ({dropped:,} lost entirely), "
          f"{exposed_days / 365.25:,.1f} exposed person-years before splitting")
    return df


# ---------------------------------------------------------------------------
# 2. Rebuilding the interval-level analysis frame
# ---------------------------------------------------------------------------
def build_lagged_frame(presc_lagged, doac_exposure, terms, term_col,
                       timeline_fn, comorbidities, binary_cols, continuous_cols,
                       ppi_pats, backfill_age=None, label=""):
    """Rebuild the interval-level Cox frame from lagged episodes.

    Mirrors, exactly, the sequence the notebook uses for the primary analysis:
    build_exposure_timeline -> attach study dates / event date / event flag ->
    add the patients with no exposure at all -> exposure indicator columns ->
    baseline covariates (NaNs preserved for MICE) -> PPI_baseline -> drop
    zero-length intervals.

    `doac_exposure` must be the outcome-specific cohort frame (patid,
    study_start_date, study_end_date, event_date).
    """
    t0 = time.time()
    cohort = doac_exposure.drop_duplicates(subset=["patid"]).copy()
    for c in ("study_start_date", "study_end_date", "event_date"):
        cohort[c] = pd.to_datetime(cohort[c], errors="coerce")
    cohort_idx = cohort.set_index("patid")

    # --- split follow-up at every exposure start/stop ----------------------
    if len(presc_lagged):
        tl = timeline_fn(presc_lagged)
    else:  # every episode was destroyed by the lag
        tl = pd.DataFrame(columns=["patid", "start_date", "end_date",
                                   "concomitant_drugs_ls"])

    if len(tl):
        tl["study_start_date"] = tl["patid"].map(cohort_idx["study_start_date"])
        tl["study_end_date"] = tl["patid"].map(cohort_idx["study_end_date"])
        tl["event_date"] = tl["patid"].map(cohort_idx["event_date"])
        tl["event_flag"] = np.where(
            (tl["event_date"] >= tl["start_date"]) &
            (tl["event_date"] <= tl["end_date"]), 1, 0)

    # --- patients contributing only unexposed time -------------------------
    have = set(tl["patid"].unique()) if len(tl) else set()
    rest = cohort[~cohort["patid"].isin(have)].copy()
    rest["concomitant_drugs_ls"] = [[] for _ in range(len(rest))]
    rest["start_date"] = rest["study_start_date"]
    rest["end_date"] = rest["study_end_date"]
    rest["event_flag"] = np.where(
        (rest["event_date"] >= rest["start_date"]) &
        (rest["event_date"] <= rest["end_date"]), 1, 0)

    df = pd.concat([tl, rest], ignore_index=True)
    df = df[df["start_date"] < df["end_date"]].reset_index(drop=True)

    # --- exposure indicators ----------------------------------------------
    for t in terms:
        df[t] = df["concomitant_drugs_ls"].apply(
            lambda v, t=t: 1 if isinstance(v, (list, tuple, set)) and t in v else 0)

    # --- baseline covariates (NaNs deliberately preserved for MICE) --------
    com = comorbidities.drop_duplicates(subset=["patid"]).set_index("patid")
    for c in binary_cols:
        if c in com.columns:
            df[c] = df["patid"].map(com[c]).fillna(0).astype(int)
    for c in continuous_cols:
        if c in com.columns:
            df[c] = pd.to_numeric(df["patid"].map(com[c]), errors="coerce")

    df["PPI_baseline"] = df["patid"].isin(ppi_pats).astype(int)

    if backfill_age is not None:
        df = backfill_age(df, doac_exposure, label)

    n_ev = int(df.loc[df["event_flag"] == 1, "patid"].nunique())
    print(f"[{label}] frame rebuilt in {time.time() - t0:,.0f}s: {len(df):,} intervals, "
          f"{df['patid'].nunique():,} patients, {n_ev:,} events")
    return df


# ---------------------------------------------------------------------------
# 3. Exposure retention diagnostics
# ---------------------------------------------------------------------------
def events_and_py(df, mask):
    """Events (unique patients) and person-years over a boolean row mask.

    total_seconds rather than .dt.days: .days truncates towards zero and
    silently shortens person-time across millions of interval rows.
    """
    sub = df.loc[mask]
    if len(sub) == 0:
        return 0, 0.0, np.nan
    days = ((sub["end_date"] - sub["start_date"]).dt.total_seconds() / 86400.0).clip(lower=0)
    py = float(days.sum()) / 365.25
    ev = int(sub.loc[sub["event_flag"] == 1, "patid"].nunique())
    return ev, py, (ev / py * 100 if py > 0 else np.nan)


def exposure_retention(frames_by_lag, terms, outcome, min_events=10):
    """Exposed events / person-years retained at each lag.

    This table is not optional. Report it beside every lagged hazard ratio so a
    reader can see when an estimate has become uninformative rather than
    reassuring. `min_events` flags rows that fall below the study's own
    pre-specified 10-event threshold.
    """
    base = frames_by_lag.get(0)
    rows = []
    for lag, df in sorted(frames_by_lag.items()):
        for t in terms:
            ev, py, rate = events_and_py(df, (df[t] == 1).to_numpy(dtype=bool))
            row = {"outcome": outcome, "lag_days": lag, "term": t,
                   "exposed_events": ev, "exposed_person_years": round(py, 1),
                   "rate_per_100py": round(rate, 2) if np.isfinite(rate) else np.nan,
                   "patients_exposed": int(df.loc[df[t] == 1, "patid"].nunique())}
            if base is not None and lag != 0:
                ev0, py0, _ = events_and_py(base, (base[t] == 1).to_numpy(dtype=bool))
                row["pct_events_retained"] = round(100 * ev / ev0, 1) if ev0 else np.nan
                row["pct_person_years_retained"] = round(100 * py / py0, 1) if py0 else np.nan
            row["below_event_threshold"] = ev < min_events
            rows.append(row)
    out = pd.DataFrame(rows)
    lost = out[(out["lag_days"] != 0) & out["below_event_threshold"]]
    if len(lost):
        print(f"[{outcome}] WARNING - below {min_events} exposed events after lagging: "
              + ", ".join(f"{r.term} (lag {r.lag_days}d: {r.exposed_events})"
                          for r in lost.itertuples()))
    return out


# ---------------------------------------------------------------------------
# 4. Driver
# ---------------------------------------------------------------------------
def run_lag_sensitivity(presc_merged, doac_exposure, imputations, *, terms, term_col,
                        ppi_adjusted, outcome, model_fns, timeline_fn,
                        comorbidities, binary_cols, continuous_cols, ppi_pats,
                        lags=None, mode="delay_onset", backfill_age=None,
                        include_crude=True, keep_frames=False,
                        results_dir="results", file_stem="lag_sensitivity"):
    """Run the full pipeline once per lag for ONE outcome cohort.

    `model_fns` is a dict of the notebook's own fitting / tabulating functions,
    so the modelling code is never duplicated here:

        model_fns = dict(
            fit_mi      = fit_adjusted_models_mi,     # (df, imputations, label=)
            summary     = build_split_ppi_summary,    # (fits, df, terms, ppi_adjusted)
            attach_diag = attach_mi_diagnostics,      # (summary, pooled, ppi_adjusted)
            crude       = crude_hr_table,             # (df, terms, label=)
        )

    Returns (summaries_by_lag, retention_table, frames_by_lag).
    `frames_by_lag` is empty unless keep_frames=True - the interval frames are
    large and are dropped as soon as the models are fitted.
    """
    import os
    lags = LAG_DAYS if lags is None else list(lags)
    os.makedirs(results_dir, exist_ok=True)

    summaries, frames, retention_frames = {}, {}, {}
    for lag in lags:
        tag = f"{outcome} lag {lag}d"
        presc_l = lag_episodes(presc_merged, lag, term_col=term_col, mode=mode, label=tag)
        df_l = build_lagged_frame(
            presc_l, doac_exposure, terms, term_col, timeline_fn,
            comorbidities, binary_cols, continuous_cols, ppi_pats,
            backfill_age=backfill_age, label=tag)

        pooled, _ = model_fns["fit_mi"](df_l, imputations, label=tag)
        summ = model_fns["summary"](pooled, df_l, terms, ppi_adjusted)
        summ = model_fns["attach_diag"](summ, pooled, ppi_adjusted)
        if include_crude:
            summ = summ.merge(model_fns["crude"](df_l, terms, label=tag),
                              on="drug_class", how="left")

        summ.insert(0, "outcome", outcome)
        summ.insert(1, "lag_days", lag)
        summ.insert(2, "lag_mode", "primary (no lag)" if lag == 0 else mode)
        summaries[lag] = summ

        # Retention needs the frames; keep only the columns it uses.
        retention_frames[lag] = df_l[["patid", "start_date", "end_date",
                                      "event_flag"] + list(terms)].copy()
        if keep_frames:
            frames[lag] = df_l
        else:
            del df_l
            gc.collect()

        path = f"{results_dir}/{file_stem}_{outcome.lower().replace(' ', '_')}_lag{lag}.xlsx"
        summ.to_excel(path, index=False)
        print(f"[{tag}] wrote {path}\n")

    retention = exposure_retention(retention_frames, terms, outcome)
    retention.to_excel(
        f"{results_dir}/{file_stem}_retention_{outcome.lower().replace(' ', '_')}.xlsx",
        index=False)
    del retention_frames
    gc.collect()
    return summaries, retention, frames


# ---------------------------------------------------------------------------
# 5. Assembling the supplementary tables
# ---------------------------------------------------------------------------
def combine_lag_summaries(summaries_by_outcome):
    """Long format: one row per outcome x lag x term. Feeds the forest plot."""
    parts = []
    for outcome, by_lag in summaries_by_outcome.items():
        for lag, summ in sorted(by_lag.items()):
            parts.append(summ)
    out = pd.concat(parts, ignore_index=True)
    return out.sort_values(["outcome", "drug_class", "lag_days"]).reset_index(drop=True)


def lag_comparison_wide(long_table, outcome=None, retention=None, hr_dp=2):
    """Wide format for the supplement: one row per drug, one block per lag.

    Columns per lag: HR (95% CI), p, and the exposed events the estimate rests
    on, so the reader can see the estimate and its fragility together.
    """
    df = long_table if outcome is None else long_table[long_table["outcome"] == outcome]
    if retention is not None:
        df = df.merge(
            retention.rename(columns={"term": "drug_class"})[
                ["outcome", "lag_days", "drug_class", "exposed_events",
                 "exposed_person_years"]],
            on=["outcome", "lag_days", "drug_class"], how="left")

    df = df.copy()
    df["HR (95% CI)"] = (df["HR"].round(hr_dp).astype(str) + " (" +
                         df["CI Lower"].round(hr_dp).astype(str) + " to " +
                         df["CI Upper"].round(hr_dp).astype(str) + ")")
    df["p"] = df["p"].map(lambda v: "<0.01" if pd.notna(v) and v < 0.01
                          else (round(v, 2) if pd.notna(v) else np.nan))

    value_cols = ["HR (95% CI)", "p"]
    if "exposed_events" in df.columns:
        value_cols += ["exposed_events"]

    wide = df.pivot_table(index=["outcome", "drug_class"], columns="lag_days",
                          values=value_cols, aggfunc="first")
    wide.columns = [f"{v} | lag {int(l)}d" if l else f"{v} | primary"
                    for v, l in wide.columns]
    _rank = {"HR (95% CI)": 0, "p": 1, "exposed_events": 2, "exposed_person_years": 3}
    order = sorted(wide.columns, key=lambda c: (
        int(c.split("lag ")[1].rstrip("d")) if "lag " in c else 0,
        _rank.get(c.split(" | ")[0], 9)))
    wide = wide[order].reset_index()

    # Order rows by the primary hazard ratio, descending.
    primary = (df[df["lag_days"] == 0]
               .set_index(["outcome", "drug_class"])["HR"].to_dict())
    wide["_sort"] = [primary.get((o, d), np.nan)
                     for o, d in zip(wide["outcome"], wide["drug_class"])]
    return (wide.sort_values(["outcome", "_sort"], ascending=[True, False])
                .drop(columns="_sort").reset_index(drop=True))
