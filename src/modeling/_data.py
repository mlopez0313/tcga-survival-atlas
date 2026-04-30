"""Shared train/test loading + leak-free preprocessing for model scripts."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from ..utils.io import load_parquet, load_pickle


def load_processed(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Load every processed parquet + the saved patient split."""
    pdir = Path(cfg["paths"]["data_processed"])
    out: Dict[str, Any] = {"y": load_parquet(pdir / "y_survival.parquet")}
    for name in ["clinical", "rna", "mutation", "cnv", "methylation", "mirna",
                 "baseline", "multimodal"]:
        path = pdir / f"X_{name}.parquet"
        if path.exists():
            out[f"X_{name}"] = load_parquet(path)
    out["split"] = load_pickle(pdir / "split.pkl")
    return out


def select_feature_set(processed: Dict[str, Any], name: str) -> pd.DataFrame:
    """Resolve `feature_set` from config to a concrete X frame."""
    key = f"X_{name}"
    if key not in processed:
        raise KeyError(
            f"Feature set '{name}' not available. Found: "
            f"{sorted(k for k in processed if k.startswith('X_'))}"
        )
    return processed[key]


def split_xy(
    X: pd.DataFrame,
    y: pd.DataFrame,
    split: Dict[str, Any],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Apply the persisted patient-id split to (X, y)."""
    train_ids = [p for p in split["train"] if p in X.index]
    test_ids = [p for p in split["test"] if p in X.index]
    return (
        X.loc[train_ids].copy(),
        y.loc[train_ids].copy(),
        X.loc[test_ids].copy(),
        y.loc[test_ids].copy(),
    )


def fit_scaler_train(
    X_train: pd.DataFrame, X_test: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, StandardScaler]:
    """Fit a StandardScaler on train ONLY, transform both."""
    scaler = StandardScaler()
    X_train_s = pd.DataFrame(
        scaler.fit_transform(X_train.values),
        index=X_train.index, columns=X_train.columns,
    ).astype(np.float32)
    X_test_s = pd.DataFrame(
        scaler.transform(X_test.values),
        index=X_test.index, columns=X_test.columns,
    ).astype(np.float32)
    return X_train_s, X_test_s, scaler
