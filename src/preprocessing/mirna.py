"""miRNA expression preprocessing.

Same pattern as RNA-seq but smaller feature set (~2k miRNAs).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd

from ..utils.logging import get_logger
from ._common import keep_primary_tumor, read_tsv, tcga_sample_to_patient, variance_filter_top_k


def preprocess_mirna(mirna_path: Path, cfg: Dict[str, Any]) -> pd.DataFrame:
    log = get_logger("preprocess_mirna", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["preprocessing"]["mirna"]

    df = read_tsv(mirna_path, index_col=0)
    log.info(f"Loaded miRNA: {df.shape}  (miRNAs x samples)")

    read_cols = [c for c in df.columns if isinstance(c, str) and c.startswith("read_count_")]
    rpm_cols = [c for c in df.columns if isinstance(c, str) and c.startswith("reads_per_million_miRNA_mapped_")]
    if read_cols or rpm_cols:
        cols = rpm_cols if rpm_cols else read_cols
        log.info(f"Detected GDC miRNA table; using {len(cols)} expression columns")
        df = df[cols]
        prefix = "reads_per_million_miRNA_mapped_" if rpm_cols else "read_count_"
        df.columns = [c[len(prefix):] if c.startswith(prefix) else c for c in cols]

    keep = keep_primary_tumor(df.columns.tolist())
    if keep:
        df = df[keep]
    df = df.T
    df.index = [tcga_sample_to_patient(s) for s in df.index]
    df.index.name = "patient_id"
    df = df.groupby(level=0).mean(numeric_only=True)
    df = df.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all")

    if spec.get("log_transform", True):
        max_v = float(np.nanmax(df.values))
        if max_v > 30:
            log.info(f"Applying log2(x+1) (max={max_v:.1f})")
            df = np.log2(df.clip(lower=0) + 1.0)

    df = variance_filter_top_k(df, int(spec.get("top_variable", 200)))
    log.info(f"Final miRNA matrix: {df.shape}")
    return df
