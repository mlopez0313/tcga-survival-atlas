"""Model-agnostic permutation importance for survival models.

We measure the drop in C-index after shuffling each feature column on
the held-out set. Works for any model that exposes a callable
`risk_fn(X) -> np.ndarray`.
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np
import pandas as pd

from .metrics import concordance_index


def permutation_importance_cindex(
    risk_fn: Callable[[pd.DataFrame], np.ndarray],
    X: pd.DataFrame,
    time: np.ndarray,
    event: np.ndarray,
    n_repeats: int = 5,
    seed: int = 42,
    columns: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Return a DataFrame with one row per feature: (importance_mean, importance_std)."""
    rng = np.random.default_rng(seed)
    base = concordance_index(time, event, risk_fn(X))
    cols = columns if columns is not None else list(X.columns)
    rows = []
    for col in cols:
        deltas = []
        for _ in range(n_repeats):
            X_perm = X.copy()
            X_perm[col] = rng.permutation(X_perm[col].to_numpy())
            c = concordance_index(time, event, risk_fn(X_perm))
            deltas.append(base - c)
        rows.append({
            "feature": col,
            "importance_mean": float(np.mean(deltas)),
            "importance_std": float(np.std(deltas)),
        })
    return pd.DataFrame(rows).sort_values("importance_mean", ascending=False).reset_index(drop=True)
