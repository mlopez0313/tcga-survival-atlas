"""Clinical metadata preprocessing.

Outputs a patient-indexed table with at minimum:
    OS_time (numeric, days), OS_event (0/1), and any of
    age, sex, stage, smoking that the Xena dump exposes.

Categorical fields are one-hot encoded later inside the merge step (so
that train-only fitting stays consistent across modalities).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd

from ..utils.logging import get_logger
from ._common import first_existing, read_tsv, tcga_sample_to_patient

# Xena harmonized survival columns are usually OS.time / OS
_SURV_TIME_CANDIDATES = ["OS.time", "OS_time", "os_time", "days_to_death", "days_to_last_follow_up"]
_SURV_EVENT_CANDIDATES = ["OS", "OS_event", "vital_status.demographic", "vital_status"]


def _coerce_event(series: pd.Series) -> pd.Series:
    """Map a vital-status / OS column to {0, 1}. NaN-safe."""
    if pd.api.types.is_numeric_dtype(series):
        return (series.astype(float) > 0).astype("Int64")
    norm = series.astype(str).str.strip().str.lower()
    mapping = {"dead": 1, "deceased": 1, "1": 1, "1.0": 1, "true": 1,
               "alive": 0, "living": 0, "0": 0, "0.0": 0, "false": 0}
    return norm.map(mapping).astype("Int64")


def _parse_stage(value: Any) -> Optional[str]:
    """Collapse free-text stage to {I, II, III, IV, NA}."""
    if pd.isna(value):
        return np.nan
    s = str(value).lower()
    if "iv" in s:
        return "IV"
    if "iii" in s:
        return "III"
    if "ii" in s:
        return "II"
    if re.search(r"\bi\b|stage i", s):
        return "I"
    return np.nan


def _load_survival(survival_path: Optional[Path], log) -> Optional[pd.DataFrame]:
    if survival_path is None or not Path(survival_path).exists():
        log.warning("No dedicated survival.tsv found; trying clinical file for OS columns.")
        return None
    df = read_tsv(survival_path)
    log.info(f"Loaded survival table: {df.shape}")
    return df


def _attach_survival(
    clinical: pd.DataFrame, survival: Optional[pd.DataFrame], log
) -> pd.DataFrame:
    """Add OS_time and OS_event columns. Uses survival.tsv when available."""
    df = clinical.copy()

    if survival is not None:
        time_col = first_existing(survival.columns, _SURV_TIME_CANDIDATES)
        event_col = first_existing(survival.columns, _SURV_EVENT_CANDIDATES)
        sample_col = first_existing(survival.columns, ["sample", "_PATIENT", "submitter_id.samples"])
        if not (time_col and event_col and sample_col):
            log.warning(
                f"Survival table missing required columns "
                f"(time={time_col}, event={event_col}, id={sample_col}); "
                "falling back to clinical."
            )
            survival = None
        else:
            surv = survival[[sample_col, time_col, event_col]].copy()
            surv["patient_id"] = surv[sample_col].map(tcga_sample_to_patient)
            surv = (
                surv.dropna(subset=["patient_id"])
                .groupby("patient_id", as_index=False)
                .agg({time_col: "max", event_col: "max"})
            )
            surv = surv.rename(columns={time_col: "OS_time", event_col: "OS_event"})
            surv["OS_event"] = _coerce_event(surv["OS_event"])
            surv["OS_time"] = pd.to_numeric(surv["OS_time"], errors="coerce")
            df = df.merge(surv[["patient_id", "OS_time", "OS_event"]], on="patient_id", how="left")
            return df

    # Fallback: derive from clinical
    event_col = first_existing(df.columns, _SURV_EVENT_CANDIDATES)
    death_col = first_existing(df.columns, ["days_to_death"])
    follow_col = first_existing(df.columns, ["days_to_last_follow_up"])
    time_col = first_existing(df.columns, ["OS.time", "OS_time", "os_time"])

    if event_col is not None and (death_col is not None or follow_col is not None):
        event = _coerce_event(df[event_col])
        death = pd.to_numeric(df[death_col], errors="coerce") if death_col is not None else pd.Series(np.nan, index=df.index)
        follow = pd.to_numeric(df[follow_col], errors="coerce") if follow_col is not None else pd.Series(np.nan, index=df.index)
        df["OS_event"] = event
        df["OS_time"] = np.where(event.fillna(0).astype(int) == 1, death, follow)
        keep_cols = [c for c in ["patient_id", "age", "sex", "stage", "smoking", "OS_time", "OS_event"] if c in df.columns]
        return df[keep_cols]

    if event_col is None or time_col is None:
        log.error("Could not locate OS_time / OS_event in clinical file. Survival columns will be NaN.")
        df["OS_time"] = np.nan
        df["OS_event"] = pd.Series(pd.array([pd.NA] * len(df), dtype="Int64"))
        return df

    df["OS_event"] = _coerce_event(df[event_col])
    df["OS_time"] = pd.to_numeric(df[time_col], errors="coerce")
    keep_cols = [c for c in ["patient_id", "age", "sex", "stage", "smoking", "OS_time", "OS_event"] if c in df.columns]
    return df[keep_cols]


def preprocess_clinical(
    clinical_path: Path,
    survival_path: Optional[Path],
    cfg: Dict[str, Any],
) -> pd.DataFrame:
    """Load + clean clinical and (separate) survival files into one table."""
    log = get_logger("preprocess_clinical", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["preprocessing"]["clinical"]["columns"]

    clin = read_tsv(clinical_path)
    log.info(f"Loaded clinical: {clin.shape}")

    pid_col = first_existing(clin.columns, spec["patient_id"] + ["submitter_id", "bcr_patient_barcode"])
    if pid_col is None:
        raise ValueError(f"Patient-ID column not found. Tried {spec['patient_id']}")
    clin["patient_id"] = clin[pid_col].map(tcga_sample_to_patient)

    age_col = first_existing(clin.columns, spec["age"])
    sex_col = first_existing(clin.columns, spec["sex"])
    stage_col = first_existing(clin.columns, spec["stage"])
    smoke_col = first_existing(clin.columns, spec.get("smoking", []))

    out_cols = {"patient_id": clin["patient_id"]}
    if age_col is not None:
        age_vals = pd.to_numeric(clin[age_col], errors="coerce")
        # Some Xena dumps store age as negative days-from-birth
        if age_vals.min(skipna=True) is not np.nan and age_vals.dropna().lt(0).mean() > 0.5:
            age_vals = (-age_vals / 365.25).round(1)
        out_cols["age"] = age_vals
    if sex_col is not None:
        out_cols["sex"] = clin[sex_col].astype(str).str.lower().replace({"nan": np.nan})
    if stage_col is not None:
        out_cols["stage"] = clin[stage_col].map(_parse_stage)
    if smoke_col is not None:
        out_cols["smoking"] = pd.to_numeric(clin[smoke_col], errors="coerce")

    df = pd.DataFrame(out_cols).drop_duplicates(subset="patient_id")
    surv_source_cols = [c for c in [pid_col, "patient_id", "days_to_death", "days_to_last_follow_up", "vital_status", "vital_status.demographic", "OS", "OS_event", "OS.time", "OS_time", "os_time"] if c in clin.columns or c == "patient_id"]
    clin_surv = clin[surv_source_cols].copy()
    if "patient_id" not in clin_surv.columns:
        clin_surv["patient_id"] = clin[pid_col].map(tcga_sample_to_patient)
    clin_surv = clin_surv.drop_duplicates(subset="patient_id")
    df = df.merge(clin_surv, on="patient_id", how="left")
    df = _attach_survival(df, _load_survival(survival_path, log), log)

    # Filter rows
    if cfg["preprocessing"]["clinical"].get("drop_missing_target", True):
        before = len(df)
        df = df.dropna(subset=["OS_time", "OS_event"])
        df = df[df["OS_time"] > 0]
        df["OS_event"] = df["OS_event"].astype(int)
        log.info(f"Filtered missing/zero survival: {before} -> {len(df)} patients")

    df = df.set_index("patient_id").sort_index()
    log.info(f"Final clinical: {df.shape}; events={int(df['OS_event'].sum())}")
    return df
