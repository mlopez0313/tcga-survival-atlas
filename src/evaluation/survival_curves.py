"""Kaplan-Meier curves + log-rank test for any model's risk scores."""
from __future__ import annotations

from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")  # headless-safe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from .metrics import risk_to_groups


def plot_km_by_risk(
    time: np.ndarray,
    event: np.ndarray,
    risk: np.ndarray,
    model_name: str,
    out_path: str | Path,
    quantile: float = 0.5,
) -> Dict[str, float]:
    """Stratify patients by median risk; plot KM curves; run log-rank test."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    groups = risk_to_groups(risk, q=quantile)
    df = pd.DataFrame({"time": time, "event": event, "group": groups})

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    kmf = KaplanMeierFitter()
    for g, label in [(0, "Low risk"), (1, "High risk")]:
        sub = df[df.group == g]
        if len(sub) == 0:
            continue
        kmf.fit(sub["time"], event_observed=sub["event"], label=f"{label} (n={len(sub)})")
        kmf.plot_survival_function(ax=ax, ci_show=True)

    lr = logrank_test(
        df.loc[df.group == 0, "time"], df.loc[df.group == 1, "time"],
        event_observed_A=df.loc[df.group == 0, "event"],
        event_observed_B=df.loc[df.group == 1, "event"],
    )
    ax.set_title(f"{model_name}  |  log-rank p = {lr.p_value:.2e}")
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival probability")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)

    return {
        "logrank_p": float(lr.p_value),
        "logrank_stat": float(lr.test_statistic),
        "n_low": int((groups == 0).sum()),
        "n_high": int((groups == 1).sum()),
    }
