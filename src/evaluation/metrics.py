"""Survival metric helpers used by every model.

Keeps the contract uniform: every model returns a 1-D numpy array of
patient-level risk scores, and every evaluator expects (time, event)
arrays + risk -> a C-index.
"""
from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd
from sksurv.metrics import concordance_index_censored


def surv_y_to_structured(y: pd.DataFrame) -> np.ndarray:
    """Convert a 2-col DataFrame [OS_event, OS_time] to scikit-survival's
    structured array `dtype=[(event, '?'), (time, '<f8')]`."""
    event = y["OS_event"].astype(bool).to_numpy()
    time = y["OS_time"].astype(float).to_numpy()
    return np.array(list(zip(event, time)), dtype=[("event", "?"), ("time", "<f8")])


def concordance_index(time: np.ndarray, event: np.ndarray, risk: np.ndarray) -> float:
    """Harrell's C-index. Higher = better. Random = 0.5."""
    event = np.asarray(event).astype(bool)
    time = np.asarray(time).astype(float)
    risk = np.asarray(risk).astype(float)
    c, *_ = concordance_index_censored(event, time, risk)
    return float(c)


def risk_to_groups(risk: np.ndarray, q: float = 0.5) -> np.ndarray:
    """Split risk into low (0) / high (1) at the q-th quantile (median by default)."""
    thr = np.quantile(risk, q)
    return (np.asarray(risk) >= thr).astype(int)


def summarize_metrics(
    name: str,
    y_train: pd.DataFrame,
    risk_train: np.ndarray,
    y_test: pd.DataFrame,
    risk_test: np.ndarray,
    extra: Dict[str, float] | None = None,
) -> Dict[str, float]:
    """Build a one-row metric dict suitable for JSON / DataFrame logging."""
    out = {
        "model": name,
        "n_train": int(len(y_train)),
        "n_test": int(len(y_test)),
        "events_train": int(y_train["OS_event"].sum()),
        "events_test": int(y_test["OS_event"].sum()),
        "c_index_train": concordance_index(
            y_train["OS_time"].to_numpy(),
            y_train["OS_event"].to_numpy(),
            risk_train,
        ),
        "c_index_test": concordance_index(
            y_test["OS_time"].to_numpy(),
            y_test["OS_event"].to_numpy(),
            risk_test,
        ),
    }
    if extra:
        out.update(extra)
    return out
