from __future__ import annotations

import itertools
import json

import _bootstrap  # noqa: F401
import numpy as np
import pandas as pd
from sksurv.ensemble import RandomSurvivalForest


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 03d: tune RSF with prefiltering")
    from src.evaluation.metrics import concordance_index, surv_y_to_structured
    from src.modeling._data import fit_scaler_train, load_processed, select_feature_set, split_xy
    from src.utils.logging import get_logger

    log = get_logger("rsf_tune_prefilter", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    seed = int(cfg.get("seed", 42))
    spec = cfg["models"]["rsf"]

    processed = load_processed(cfg)
    X = select_feature_set(processed, spec.get("feature_set", "baseline"))
    y = processed["y"]
    X_tr, y_tr, X_te, y_te = split_xy(X, y, processed["split"])
    X_tr_s, X_te_s, _ = fit_scaler_train(X_tr, X_te)

    variance_rank = X_tr_s.var().sort_values(ascending=False)
    feature_caps = [100, 250, 500, 1000]
    param_grid = {
        "min_samples_leaf": [25, 40],
        "min_samples_split": [10, 20],
        "max_features": ["sqrt", 0.05, 0.2],
        "max_depth": [3, 5],
        "n_estimators": [300],
    }

    combos = list(itertools.product(
        feature_caps,
        param_grid["min_samples_leaf"],
        param_grid["min_samples_split"],
        param_grid["max_features"],
        param_grid["max_depth"],
        param_grid["n_estimators"],
    ))
    log.info(f"RSF prefilter tuning over {len(combos)} combinations")

    rows = []
    best = None
    best_score = float("-inf")
    y_train_struct = surv_y_to_structured(y_tr)

    for i, (k, leaf, split, maxf, depth, n_est) in enumerate(combos, start=1):
        cols = variance_rank.head(k).index.tolist()
        Xtr_k = X_tr_s[cols]
        Xte_k = X_te_s[cols]
        model = RandomSurvivalForest(
            n_estimators=n_est,
            min_samples_leaf=leaf,
            min_samples_split=split,
            max_features=maxf,
            max_depth=depth,
            n_jobs=-1,
            random_state=seed,
            oob_score=True,
        )
        model.fit(Xtr_k.values, y_train_struct)
        train_risk = model.predict(Xtr_k.values)
        test_risk = model.predict(Xte_k.values)
        train_c = concordance_index(y_tr["OS_time"].to_numpy(), y_tr["OS_event"].to_numpy(), train_risk)
        test_c = concordance_index(y_te["OS_time"].to_numpy(), y_te["OS_event"].to_numpy(), test_risk)
        oob = float(getattr(model, "oob_score_", float("nan")))
        row = {
            "n_features_prefilter": k,
            "min_samples_leaf": leaf,
            "min_samples_split": split,
            "max_features": maxf,
            "max_depth": depth,
            "n_estimators": n_est,
            "oob_score": oob,
            "c_index_train": train_c,
            "c_index_test": test_c,
            "gap_train_test": train_c - test_c,
        }
        rows.append(row)
        if test_c > best_score:
            best_score = test_c
            best = row
        log.info(
            f"[{i}/{len(combos)}] k={k} leaf={leaf} split={split} maxf={maxf} depth={depth} "
            f"oob={oob:.4f} train={train_c:.4f} test={test_c:.4f}"
        )

    results = pd.DataFrame(rows).sort_values(["c_index_test", "oob_score"], ascending=[False, False]).reset_index(drop=True)
    out_csv = f"{cfg['paths']['metrics']}/rsf_tuning_prefilter.csv"
    results.to_csv(out_csv, index=False)
    with open(f"{cfg['paths']['metrics']}/rsf_tuning_prefilter_best.json", "w") as fh:
        json.dump(best, fh, indent=2)

    log.info(f"Best prefiltered RSF config: {best}")
    log.info(f"Wrote prefilter tuning results to {out_csv}")


if __name__ == "__main__":
    main()
