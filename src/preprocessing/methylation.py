"""DNA methylation preprocessing (Illumina 450K beta values).

Methylation matrices are huge; we drop NaN-heavy probes early, then keep
the top-variable probes. Probe-to-gene aggregation is left as a future
enhancement (would require an annotation file).
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import pandas as pd

from ..utils.logging import get_logger
from ._common import keep_primary_tumor, read_tsv, tcga_sample_to_patient, variance_filter_top_k


def preprocess_methylation(meth_path: Path, cfg: Dict[str, Any]) -> pd.DataFrame:
    log = get_logger("preprocess_methylation", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["preprocessing"]["methylation"]

    df = read_tsv(meth_path, index_col=0)
    log.info(f"Loaded methylation: {df.shape}  (probes x samples)")

    keep = keep_primary_tumor(df.columns.tolist())
    if keep:
        df = df[keep]

    # Drop probes with too many NaNs (computed across samples = rows post-T = no, here rows are probes)
    max_missing = float(spec.get("max_missing_frac", 0.2))
    nan_frac = df.isna().mean(axis=1)
    keep_probes = nan_frac <= max_missing
    log.info(f"Probe NaN filter (<= {max_missing*100:.0f}% missing): "
             f"kept {int(keep_probes.sum())}/{len(keep_probes)}")
    df = df.loc[keep_probes]

    df = df.T
    df.index = [tcga_sample_to_patient(s) for s in df.index]
    df.index.name = "patient_id"
    df = df.groupby(level=0).mean()
    df = df.apply(pd.to_numeric, errors="coerce")

    df = variance_filter_top_k(df, int(spec.get("top_variable_probes", 5000)))
    df = df.fillna(df.median(axis=0))  # remaining NaNs -> probe median
    log.info(f"Final methylation matrix: {df.shape}")
    return df
