"""Copy-number alteration preprocessing.

Expects a Xena GISTIC matrix: (genes x samples) of integer focal calls
{-2,-1,0,1,2} or continuous segment means. We collapse to patients,
optionally binarize, and select top-variable features.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd

from ..utils.logging import get_logger
from ._common import keep_primary_tumor, read_tsv, tcga_sample_to_patient, variance_filter_top_k


def preprocess_cnv(cnv_path: Path, cfg: Dict[str, Any]) -> pd.DataFrame:
    log = get_logger("preprocess_cnv", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["preprocessing"]["cnv"]

    df = read_tsv(cnv_path)
    log.info(f"Loaded CNV: {df.shape}")

    if {"gene_id", "gene_name"}.issubset(df.columns):
        gene_col = "gene_name"
        value_cols = [c for c in df.columns if c.endswith("_copy_number")]
        log.info(f"Detected GDC gene-level CNV table with {len(value_cols)} sample columns")
        mat = df[[gene_col] + value_cols].copy()
        mat = mat.groupby(gene_col, as_index=True)[value_cols].mean()
        keep = keep_primary_tumor(value_cols)
        if keep:
            mat = mat[keep]
        mat.columns = [c.replace("_copy_number", "") for c in mat.columns]
        df = mat
    else:
        df = read_tsv(cnv_path, index_col=0)
        keep = keep_primary_tumor(df.columns.tolist())
        if keep:
            df = df[keep]

    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.T
    df.index = [tcga_sample_to_patient(s) for s in df.index]
    df.index.name = "patient_id"
    df = df.groupby(level=0).mean(numeric_only=True)
    df = df.dropna(axis=1, how="all").fillna(0.0)

    mode = spec.get("cnv_mode", "continuous")
    if mode == "binary":
        thr = float(spec.get("binary_threshold", 0.3))
        df = (df.abs() >= thr).astype("int8")
        log.info(f"Binarized CNV at |x| >= {thr}; mean events/patient = {df.sum(axis=1).mean():.1f}")
    else:
        log.info("Using continuous CNV values")

    df = variance_filter_top_k(df, int(spec.get("top_variable_features", 2000)))
    log.info(f"Final CNV matrix: {df.shape}")
    return df
