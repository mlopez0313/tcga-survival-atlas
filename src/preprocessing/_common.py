"""Helpers shared across modality preprocessors."""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd


def read_tsv(path: str | Path, **kwargs) -> pd.DataFrame:
    """Read a (possibly gzipped) TSV with sane defaults."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t", low_memory=False, **kwargs)


def tcga_sample_to_patient(sample: str) -> str:
    """`TCGA-XX-YYYY-01A-...` -> `TCGA-XX-YYYY` (patient barcode).

    Handles bare patient ids too.
    """
    if not isinstance(sample, str):
        return sample
    parts = sample.split("-")
    if len(parts) < 3:
        return sample
    return "-".join(parts[:3])


def keep_primary_tumor(samples: Iterable[str]) -> list[str]:
    """Filter Xena sample IDs down to primary-tumor aliquots (codes 01-09).

    Sample-type code is the 4th block of the barcode (positions 0-1 numeric).
    Code 01-09 = tumor; 10-19 = normal; 20+ = control. We keep tumor.
    """
    out = []
    for s in samples:
        if not isinstance(s, str):
            continue
        parts = s.split("-")
        if len(parts) < 4:
            out.append(s)
            continue
        code = parts[3][:2]
        if code.isdigit() and 1 <= int(code) <= 9:
            out.append(s)
    return out


def first_existing(df_columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    """Return the first column in `candidates` that exists in `df_columns`."""
    cols = set(df_columns)
    for c in candidates:
        if c in cols:
            return c
    return None


def variance_filter_top_k(df: pd.DataFrame, k: int) -> pd.DataFrame:
    """Keep the `k` columns with highest variance."""
    if df.shape[1] <= k:
        return df
    var = df.var(axis=0, skipna=True).sort_values(ascending=False)
    return df.loc[:, var.head(k).index]
