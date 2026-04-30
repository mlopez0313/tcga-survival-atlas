"""Filesystem / serialization helpers.

Small wrappers used by every script so paths and formats stay consistent.
"""
from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any, Dict

import pandas as pd
import yaml


def load_config(path: str | Path) -> Dict[str, Any]:
    """Load a YAML config file into a plain dict."""
    with open(path, "r") as fh:
        cfg = yaml.safe_load(fh)
    if cfg is None:
        raise ValueError(f"Config at {path} is empty.")
    return cfg


def ensure_dirs(cfg: Dict[str, Any]) -> Dict[str, Path]:
    """Materialize every path under cfg['paths'] and return resolved Paths."""
    paths = {k: Path(v) for k, v in cfg["paths"].items()}
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    return paths


def save_pickle(obj: Any, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        pickle.dump(obj, fh)


def load_pickle(path: str | Path) -> Any:
    with open(path, "rb") as fh:
        return pickle.load(fh)


def save_parquet(df: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=True)


def load_parquet(path: str | Path) -> pd.DataFrame:
    return pd.read_parquet(path)
