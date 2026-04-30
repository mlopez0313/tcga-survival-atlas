"""Cox elastic-net survival model (scikit-survival)."""
from __future__ import annotations

import json
import math
import warnings
from pathlib import Path
from typing import Any, Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis

from ..evaluation.metrics import concordance_index, summarize_metrics, surv_y_to_structured
from ..evaluation.survival_curves import plot_km_by_risk
from ..utils.io import save_pickle
from ..utils.logging import get_logger
from ._data import fit_scaler_train, load_processed, select_feature_set, split_xy


def _cv_pick_alpha(X: np.ndarray, y_struct: np.ndarray, l1_ratio: float, n_alphas: int,
                   cv: int, log) -> float:
    """Two-stage: fit a path for the alpha grid, then K-fold pick the best."""
    base = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, n_alphas=n_alphas, alpha_min_ratio=0.01)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ConvergenceWarning)
        base.fit(X, y_struct)
    alphas = base.alphas_
    log.info(f"Tuning over {len(alphas)} alphas via {cv}-fold CV (l1_ratio={l1_ratio})")

    kf = KFold(n_splits=cv, shuffle=True, random_state=0)
    sums = {a: 0.0 for a in alphas}
    counts = {a: 0 for a in alphas}
    for fold, (tr, va) in enumerate(kf.split(X)):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=ConvergenceWarning)
                m = CoxnetSurvivalAnalysis(l1_ratio=l1_ratio, alphas=list(alphas),
                                           fit_baseline_model=False)
                m.fit(X[tr], y_struct[tr])
        except Exception as exc:  # noqa: BLE001
            log.warning(f"fold {fold} failed: {exc}")
            continue
        for a in m.alphas_:
            try:
                risk = m.predict(X[va], alpha=a)
                c = concordance_index(y_struct["time"][va], y_struct["event"][va], risk)
            except Exception:
                continue
            # Match back to the closest alpha in the master grid
            key = min(sums, key=lambda x: abs(x - a))
            sums[key] += c
            counts[key] += 1
    means = {a: (sums[a] / counts[a]) for a in alphas if counts[a] > 0}
    if not means:
        log.warning("CV produced no usable scores; falling back to median alpha")
        return float(alphas[len(alphas) // 2])
    best_alpha = max(means, key=means.get)
    log.info(f"Best alpha = {best_alpha:.4g} (CV C-index = {means[best_alpha]:.4f})")
    return float(best_alpha)


def _fit_with_alpha_backoff(X: np.ndarray, y_struct: np.ndarray, *, l1_ratio: float,
                            start_alpha: float, log) -> tuple[CoxnetSurvivalAnalysis, float]:
    alpha = float(start_alpha)
    multipliers = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    last_exc = None
    for mult in multipliers:
        trial_alpha = alpha * mult
        model = CoxnetSurvivalAnalysis(
            l1_ratio=l1_ratio,
            alphas=[trial_alpha],
            fit_baseline_model=True,
        )
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=ConvergenceWarning)
                model.fit(X, y_struct)
            if mult > 1.0:
                log.warning(f"Cox fit required alpha backoff: {alpha:.4g} -> {trial_alpha:.4g}")
            return model, float(trial_alpha)
        except ArithmeticError as exc:
            last_exc = exc
            log.warning(f"Cox fit failed at alpha={trial_alpha:.4g}: {exc}")
            continue
    raise last_exc if last_exc is not None else RuntimeError("Cox fit failed for unknown reason")


def train_cox_elastic_net(cfg: Dict[str, Any]) -> Dict[str, Any]:
    log = get_logger("cox_elastic_net", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    spec = cfg["models"]["cox_elastic_net"]

    processed = load_processed(cfg)
    X = select_feature_set(processed, spec.get("feature_set", "baseline"))
    y = processed["y"]
    X_tr, y_tr, X_te, y_te = split_xy(X, y, processed["split"])
    X_tr_s, X_te_s, scaler = fit_scaler_train(X_tr, X_te)

    y_train_struct = surv_y_to_structured(y_tr)
    y_test_struct = surv_y_to_structured(y_te)

    best_alpha = _cv_pick_alpha(
        X_tr_s.values.astype(np.float64),
        y_train_struct,
        l1_ratio=float(spec.get("l1_ratio", 0.5)),
        n_alphas=int(spec.get("n_alphas", 50)),
        cv=int(spec.get("cv_folds", 5)),
        log=log,
    )

    model, best_alpha = _fit_with_alpha_backoff(
        X_tr_s.values.astype(np.float64),
        y_train_struct,
        l1_ratio=float(spec.get("l1_ratio", 0.5)),
        start_alpha=best_alpha,
        log=log,
    )

    risk_train = model.predict(X_tr_s.values.astype(np.float64))
    risk_test = model.predict(X_te_s.values.astype(np.float64))

    metrics = summarize_metrics(
        "cox_elastic_net", y_tr, risk_train, y_te, risk_test,
        extra={"alpha": best_alpha, "l1_ratio": float(spec.get("l1_ratio", 0.5)),
               "n_features": int(X_tr.shape[1]),
               "feature_set": spec.get("feature_set", "baseline")},
    )
    log.info(f"Cox elastic-net   C-index train={metrics['c_index_train']:.4f}  test={metrics['c_index_test']:.4f}")

    # ------------------------------------------------------------------
    # Persist artifacts
    # ------------------------------------------------------------------
    models_dir = Path(cfg["paths"]["models"])
    fig_dir = Path(cfg["paths"]["figures"])
    metrics_dir = Path(cfg["paths"]["metrics"])
    save_pickle(
        {"model": model, "scaler": scaler, "feature_names": list(X_tr.columns),
         "alpha": best_alpha, "feature_set": spec.get("feature_set", "baseline")},
        models_dir / "cox_elastic_net.pkl",
    )

    # KM stratification on test
    km = plot_km_by_risk(
        y_te["OS_time"].to_numpy(), y_te["OS_event"].to_numpy(),
        risk_test, "Cox elastic-net (test)",
        fig_dir / "km_cox_elastic_net.png",
    )
    metrics.update({f"km_test_{k}": v for k, v in km.items()})

    # Top coefficients plot
    coef = pd.Series(model.coef_.ravel(), index=X_tr.columns)
    coef_nonzero = coef[coef != 0].sort_values()
    if len(coef_nonzero):
        top_n = min(20, len(coef_nonzero))
        top = pd.concat([coef_nonzero.head(top_n // 2), coef_nonzero.tail(top_n // 2)])
        fig, ax = plt.subplots(figsize=(7, max(4, 0.3 * len(top))))
        colors = ["#d62728" if v > 0 else "#1f77b4" for v in top.values]
        ax.barh(range(len(top)), top.values, color=colors)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top.index, fontsize=8)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_xlabel("Cox coefficient (positive = higher risk)")
        ax.set_title(f"Cox elastic-net top features (alpha={best_alpha:.3g})")
        fig.tight_layout()
        fig.savefig(fig_dir / "cox_top_coefficients.png", dpi=150)
        plt.close(fig)
        coef_nonzero.to_csv(metrics_dir / "cox_nonzero_coefficients.csv", header=["coef"])

    with open(metrics_dir / "cox_elastic_net.json", "w") as fh:
        json.dump(metrics, fh, indent=2)

    return {
        "model": model, "scaler": scaler, "metrics": metrics,
        "risk_train": risk_train, "risk_test": risk_test,
        "X_train": X_tr_s, "X_test": X_te_s, "y_train": y_tr, "y_test": y_te,
    }
