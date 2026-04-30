"""Random Survival Forest baseline (scikit-survival)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sksurv.ensemble import RandomSurvivalForest

from ..evaluation.metrics import summarize_metrics, surv_y_to_structured
from ..evaluation.survival_curves import plot_km_by_risk
from ..utils.io import save_pickle
from ..utils.logging import get_logger
from ._data import fit_scaler_train, load_processed, select_feature_set, split_xy


def _permutation_importance(model, X_te: pd.DataFrame, y_te_struct, n_repeats: int, seed: int, log, cols: list[str] | None = None) -> pd.DataFrame:
    """Permutation importance against the model's `predict` (cumulative hazard)."""
    rng = np.random.default_rng(seed)
    base = model.score(X_te.values, y_te_struct)
    rows = []
    cols = cols or list(X_te.columns)
    log.info(f"Permutation importance over {len(cols)} features ({n_repeats} repeats)...")
    for col in cols:
        deltas = []
        col_idx = X_te.columns.get_loc(col)
        Xv = X_te.values.copy()
        original = Xv[:, col_idx].copy()
        for _ in range(n_repeats):
            Xv[:, col_idx] = rng.permutation(original)
            try:
                c = model.score(Xv, y_te_struct)
            except Exception:
                c = base
            deltas.append(base - c)
        Xv[:, col_idx] = original
        rows.append({"feature": col, "importance_mean": float(np.mean(deltas)),
                     "importance_std": float(np.std(deltas))})
    return pd.DataFrame(rows).sort_values("importance_mean", ascending=False).reset_index(drop=True)


def train_rsf(cfg: Dict[str, Any]) -> Dict[str, Any]:
    log = get_logger("rsf", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["models"]["rsf"]
    seed = int(cfg.get("seed", 42))

    processed = load_processed(cfg)
    X = select_feature_set(processed, spec.get("feature_set", "baseline"))
    y = processed["y"]
    X_tr, y_tr, X_te, y_te = split_xy(X, y, processed["split"])
    X_tr_s, X_te_s, scaler = fit_scaler_train(X_tr, X_te)

    y_train_struct = surv_y_to_structured(y_tr)
    y_test_struct = surv_y_to_structured(y_te)

    model = RandomSurvivalForest(
        n_estimators=int(spec.get("n_estimators", 300)),
        min_samples_split=int(spec.get("min_samples_split", 10)),
        min_samples_leaf=int(spec.get("min_samples_leaf", 15)),
        max_features=spec.get("max_features", "sqrt"),
        n_jobs=-1,
        random_state=seed,
    )
    log.info(f"Fitting RSF on X_train: {X_tr_s.shape}")
    model.fit(X_tr_s.values, y_train_struct)

    risk_train = model.predict(X_tr_s.values)
    risk_test = model.predict(X_te_s.values)

    metrics = summarize_metrics(
        "random_survival_forest", y_tr, risk_train, y_te, risk_test,
        extra={"n_estimators": int(spec.get("n_estimators", 300)),
               "n_features": int(X_tr.shape[1]),
               "feature_set": spec.get("feature_set", "baseline")},
    )
    log.info(f"RSF   C-index train={metrics['c_index_train']:.4f}  test={metrics['c_index_test']:.4f}")

    # Persist
    models_dir = Path(cfg["paths"]["models"])
    fig_dir = Path(cfg["paths"]["figures"])
    metrics_dir = Path(cfg["paths"]["metrics"])

    save_pickle(
        {"model": model, "scaler": scaler, "feature_names": list(X_tr.columns),
         "feature_set": spec.get("feature_set", "baseline")},
        models_dir / "rsf.pkl",
    )

    km = plot_km_by_risk(
        y_te["OS_time"].to_numpy(), y_te["OS_event"].to_numpy(),
        risk_test, "RSF (test)", fig_dir / "km_rsf.png",
    )
    metrics.update({f"km_test_{k}": v for k, v in km.items()})

    # Permutation importance (cap permuted features for speed, but keep full model input width)
    cap = 200
    X_te_imp = X_te_s.copy()
    if X_te_s.shape[1] > cap:
        log.info(f"Importance capped to top-{cap} highest-variance features for speed")
        var = X_tr_s.var().sort_values(ascending=False)
        importance_cols = var.head(cap).index.tolist()
    else:
        importance_cols = list(X_te_s.columns)

    importance = _permutation_importance(model, X_te_imp, y_test_struct, n_repeats=3, seed=seed, log=log, cols=importance_cols)
    importance.to_csv(metrics_dir / "rsf_permutation_importance.csv", index=False)

    top = importance.head(20)
    fig, ax = plt.subplots(figsize=(7, max(4, 0.3 * len(top))))
    ax.barh(range(len(top))[::-1], top["importance_mean"], color="#2ca02c")
    ax.set_yticks(range(len(top))[::-1])
    ax.set_yticklabels(top["feature"], fontsize=8)
    ax.set_xlabel("Permutation importance (Δ C-index)")
    ax.set_title("RSF top features")
    fig.tight_layout()
    fig.savefig(fig_dir / "rsf_top_features.png", dpi=150)
    plt.close(fig)

    with open(metrics_dir / "rsf.json", "w") as fh:
        json.dump(metrics, fh, indent=2)

    return {
        "model": model, "scaler": scaler, "metrics": metrics,
        "risk_train": risk_train, "risk_test": risk_test,
        "X_train": X_tr_s, "X_test": X_te_s, "y_train": y_tr, "y_test": y_te,
    }
