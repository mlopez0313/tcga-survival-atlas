from __future__ import annotations

import itertools
import json

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sksurv.ensemble import RandomSurvivalForest


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 03c: tune Random Survival Forest with CV")
    from src.evaluation.metrics import concordance_index, surv_y_to_structured
    from src.modeling._data import fit_scaler_train, load_processed, select_feature_set, split_xy
    from src.utils.logging import get_logger

    log = get_logger("rsf_tune_cv", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    seed = int(cfg.get("seed", 42))
    spec = cfg["models"]["rsf"]

    processed = load_processed(cfg)
    X = select_feature_set(processed, spec.get("feature_set", "baseline"))
    y = processed["y"]
    X_tr, y_tr, X_te, y_te = split_xy(X, y, processed["split"])
    X_tr_s, X_te_s, _ = fit_scaler_train(X_tr, X_te)

    grid = {
        "min_samples_leaf": [10, 25, 40],
        "min_samples_split": [10, 20, 40],
        "max_features": ["sqrt", 0.05, 0.2],
        "max_depth": [3, 5, None],
        "n_estimators": [300],
    }

    combos = list(itertools.product(
        grid["min_samples_leaf"],
        grid["min_samples_split"],
        grid["max_features"],
        grid["max_depth"],
        grid["n_estimators"],
    ))
    log.info(f"RSF CV tuning over {len(combos)} combinations")

    cv = KFold(n_splits=5, shuffle=True, random_state=seed)
    rows = []
    best = None
    best_cv = float("-inf")

    for i, (leaf, split, maxf, depth, n_est) in enumerate(combos, start=1):
        fold_scores = []
        for tr_idx, va_idx in cv.split(X_tr_s):
            X_fold_tr = X_tr_s.iloc[tr_idx]
            X_fold_va = X_tr_s.iloc[va_idx]
            y_fold_tr = y_tr.iloc[tr_idx]
            y_fold_va = y_tr.iloc[va_idx]
            model = RandomSurvivalForest(
                n_estimators=n_est,
                min_samples_leaf=leaf,
                min_samples_split=split,
                max_features=maxf,
                max_depth=depth,
                n_jobs=-1,
                random_state=seed,
            )
            model.fit(X_fold_tr.values, surv_y_to_structured(y_fold_tr))
            risk_va = model.predict(X_fold_va.values)
            c_va = concordance_index(y_fold_va["OS_time"].to_numpy(), y_fold_va["OS_event"].to_numpy(), risk_va)
            fold_scores.append(c_va)

        cv_mean = float(np.mean(fold_scores))
        cv_std = float(np.std(fold_scores))

        final_model = RandomSurvivalForest(
            n_estimators=n_est,
            min_samples_leaf=leaf,
            min_samples_split=split,
            max_features=maxf,
            max_depth=depth,
            n_jobs=-1,
            random_state=seed,
            oob_score=True,
        )
        final_model.fit(X_tr_s.values, surv_y_to_structured(y_tr))
        train_risk = final_model.predict(X_tr_s.values)
        test_risk = final_model.predict(X_te_s.values)
        train_c = concordance_index(y_tr["OS_time"].to_numpy(), y_tr["OS_event"].to_numpy(), train_risk)
        test_c = concordance_index(y_te["OS_time"].to_numpy(), y_te["OS_event"].to_numpy(), test_risk)
        oob = float(getattr(final_model, "oob_score_", float("nan")))

        row = {
            "min_samples_leaf": leaf,
            "min_samples_split": split,
            "max_features": maxf,
            "max_depth": depth,
            "n_estimators": n_est,
            "cv_mean": cv_mean,
            "cv_std": cv_std,
            "oob_score": oob,
            "c_index_train": train_c,
            "c_index_test": test_c,
            "gap_train_test": train_c - test_c,
        }
        rows.append(row)
        if cv_mean > best_cv:
            best_cv = cv_mean
            best = row
        log.info(
            f"[{i}/{len(combos)}] leaf={leaf} split={split} maxf={maxf} depth={depth} "
            f"cv={cv_mean:.4f}±{cv_std:.4f} oob={oob:.4f} train={train_c:.4f} test={test_c:.4f}"
        )

    results = pd.DataFrame(rows).sort_values(["cv_mean", "c_index_test"], ascending=[False, False]).reset_index(drop=True)
    out_csv = f"{cfg['paths']['metrics']}/rsf_tuning_cv.csv"
    results.to_csv(out_csv, index=False)
    with open(f"{cfg['paths']['metrics']}/rsf_tuning_cv_best.json", "w") as fh:
        json.dump(best, fh, indent=2)

    log.info(f"Best CV-selected RSF config: {best}")
    log.info(f"Wrote CV tuning results to {out_csv}")


if __name__ == "__main__":
    main()
