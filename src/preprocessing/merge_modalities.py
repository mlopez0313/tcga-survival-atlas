"""Align modalities on the shared set of patient IDs and emit train/test splits.

This module ONLY:
    1. intersects patient IDs across modalities,
    2. one-hot encodes clinical categoricals,
    3. assembles `baseline` (clinical+RNA) and `multimodal` (everything),
    4. produces a single train/test split that every downstream model uses.

Imputation, scaling, and feature selection that depend on data values
are LEFT to the model-training scripts so they can be fit on train only,
keeping the pipeline leak-free.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from ..utils.io import save_parquet, save_pickle
from ..utils.logging import get_logger

# Categorical clinical fields we one-hot encode (numeric ones pass through)
_CATEGORICALS = ["sex", "stage"]


def _encode_clinical(clinical: pd.DataFrame, log) -> pd.DataFrame:
    df = clinical.copy()
    survival = df[["OS_time", "OS_event"]]
    feats = df.drop(columns=["OS_time", "OS_event"])

    cat_cols = [c for c in _CATEGORICALS if c in feats.columns]
    num_cols = [c for c in feats.columns if c not in cat_cols]

    if cat_cols:
        feats = pd.get_dummies(feats, columns=cat_cols, drop_first=False, dummy_na=False)
        # bool dummies -> int8
        for c in feats.columns:
            if feats[c].dtype == bool:
                feats[c] = feats[c].astype("int8")
    log.info(f"Clinical features after one-hot: {feats.shape[1]} cols ({cat_cols} encoded)")
    return pd.concat([feats, survival], axis=1)


def _intersect(modalities: Dict[str, pd.DataFrame]) -> list[str]:
    idx = None
    for name, df in modalities.items():
        if df is None:
            continue
        s = set(df.index)
        idx = s if idx is None else idx & s
    return sorted(idx) if idx else []


def merge_modalities(
    clinical: pd.DataFrame,
    rna: Optional[pd.DataFrame],
    mutation: Optional[pd.DataFrame],
    cnv: Optional[pd.DataFrame],
    methylation: Optional[pd.DataFrame],
    mirna: Optional[pd.DataFrame],
    cfg: Dict[str, Any],
) -> Dict[str, Any]:
    """Align modalities, emit train/test split, and persist to data/processed/."""
    log = get_logger("merge_modalities", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))

    clinical_enc = _encode_clinical(clinical, log)
    available = {
        "clinical": clinical_enc.drop(columns=["OS_time", "OS_event"]),
        "rna": rna,
        "mutation": mutation,
        "cnv": cnv,
        "methylation": methylation,
        "mirna": mirna,
    }
    available = {k: v for k, v in available.items() if v is not None and not v.empty}
    log.info(f"Modalities present: {list(available.keys())}")

    # Patient intersection across (clinical + survival) and every other modality.
    survival_index = clinical_enc.dropna(subset=["OS_time", "OS_event"]).index
    patient_ids = sorted(set(survival_index).intersection(*[set(df.index) for df in available.values()]))
    if not patient_ids:
        raise RuntimeError("No patients shared across modalities. Check sample-id parsing.")
    log.info(f"Shared patients across modalities: {len(patient_ids)}")

    # Reindex everything to the shared patient list
    aligned: Dict[str, pd.DataFrame] = {
        name: df.reindex(patient_ids) for name, df in available.items()
    }
    survival = clinical_enc.loc[patient_ids, ["OS_time", "OS_event"]].copy()
    survival["OS_event"] = survival["OS_event"].astype(int)

    # Numeric clinical NaNs are not allowed by survival models -> simple median impute
    # (training-time scalers fit on train and re-apply, but a pre-impute keeps types numeric).
    clin = aligned["clinical"]
    num_cols = clin.select_dtypes(include=[np.number, "Int64"]).columns
    for c in num_cols:
        clin[c] = pd.to_numeric(clin[c], errors="coerce")
    clin = clin.fillna(clin.median(numeric_only=True))
    # Drop any remaining all-NaN columns (e.g. if smoking missing entirely)
    clin = clin.dropna(axis=1, how="all")
    aligned["clinical"] = clin.astype(np.float32)

    # Concatenated views used by single-input models
    rna_part = aligned.get("rna")
    baseline_parts = [aligned["clinical"]]
    if rna_part is not None:
        baseline_parts.append(rna_part.add_prefix("rna__"))
    X_baseline = pd.concat(baseline_parts, axis=1).astype(np.float32)

    mm_parts = [aligned["clinical"]]
    for name, prefix in [
        ("rna", "rna__"),
        ("mutation", "mut__"),
        ("cnv", "cnv__"),
        ("methylation", "meth__"),
        ("mirna", "mir__"),
    ]:
        if name in aligned:
            mm_parts.append(aligned[name].add_prefix(prefix))
    X_multimodal = pd.concat(mm_parts, axis=1).astype(np.float32)

    log.info(f"X_baseline:   {X_baseline.shape}")
    log.info(f"X_multimodal: {X_multimodal.shape}")

    # ------------------------------------------------------------------
    # Train / test split (single split shared by every model)
    # ------------------------------------------------------------------
    seed = int(cfg.get("seed", 42))
    test_size = float(cfg["split"]["test_size"])
    stratify = survival["OS_event"] if cfg["split"].get("stratify_by_event", True) else None
    train_ids, test_ids = train_test_split(
        patient_ids, test_size=test_size, random_state=seed, stratify=stratify
    )
    log.info(
        f"Split: train={len(train_ids)} (events={int(survival.loc[train_ids,'OS_event'].sum())}) "
        f"test={len(test_ids)} (events={int(survival.loc[test_ids,'OS_event'].sum())})"
    )

    # ------------------------------------------------------------------
    # Persist everything to data/processed/
    # ------------------------------------------------------------------
    out_dir = Path(cfg["paths"]["data_processed"])
    out_dir.mkdir(parents=True, exist_ok=True)

    for name, df in aligned.items():
        save_parquet(df, out_dir / f"X_{name}.parquet")
    save_parquet(X_baseline, out_dir / "X_baseline.parquet")
    save_parquet(X_multimodal, out_dir / "X_multimodal.parquet")
    save_parquet(survival, out_dir / "y_survival.parquet")

    split = {"train": list(train_ids), "test": list(test_ids), "seed": seed}
    save_pickle(split, out_dir / "split.pkl")
    log.info(f"Wrote processed tables to {out_dir}")

    return {
        "modalities": aligned,
        "X_baseline": X_baseline,
        "X_multimodal": X_multimodal,
        "y": survival,
        "split": split,
    }
