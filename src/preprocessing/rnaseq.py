"""RNA-seq preprocessing.

Xena gene-expression matrices are typically (genes x samples) with the
gene id in the first column. We:
    1. transpose to (samples x genes),
    2. collapse multiple samples per patient (mean of tumor aliquots),
    3. low-expression filter,
    4. log2(x+1) if requested,
    5. variance top-K,
returning a (patients x genes) DataFrame indexed by patient barcode.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd

from ..utils.logging import get_logger
from ._common import keep_primary_tumor, read_tsv, tcga_sample_to_patient, variance_filter_top_k


def preprocess_rnaseq(rnaseq_path: Path, cfg: Dict[str, Any]) -> pd.DataFrame:
    log = get_logger("preprocess_rnaseq", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["preprocessing"]["rnaseq"]

    df = read_tsv(rnaseq_path, index_col=0)
    log.info(f"Loaded RNA-seq matrix: {df.shape}  (genes x samples)")

    tumor_samples = keep_primary_tumor(df.columns.tolist())
    if tumor_samples:
        df = df[tumor_samples]
    log.info(f"Kept primary-tumor samples: {df.shape}")

    df = df.T  # samples x genes
    df.index = [tcga_sample_to_patient(s) for s in df.index]
    df.index.name = "patient_id"
    df = df.groupby(level=0).mean()  # collapse multiple aliquots
    log.info(f"Collapsed to patients: {df.shape}")

    df = df.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all")

    # Detect whether values are already log-transformed: if the max is > ~30,
    # the data is raw counts/FPKM and benefits from log2(x+1).
    if spec.get("log_transform", True):
        max_v = float(np.nanmax(df.values))
        if max_v > 30:
            log.info(f"Applying log2(x+1) (max value before = {max_v:.1f})")
            df = np.log2(df.clip(lower=0) + 1.0)
        else:
            log.info(f"Skipping log transform (max={max_v:.2f} suggests already log-scale)")

    if (mn := spec.get("min_expression")) is not None:
        keep = df.mean(axis=0) >= mn
        log.info(f"Min-expression filter (>= {mn}): kept {int(keep.sum())}/{len(keep)} genes")
        df = df.loc[:, keep]

    df = variance_filter_top_k(df, int(spec.get("top_variable_genes", 2000)))
    log.info(f"Top-variable selection: {df.shape}")
    return df
