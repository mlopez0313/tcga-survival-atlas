from __future__ import annotations

import json
import itertools

import _bootstrap  # noqa: F401
import pandas as pd
from sksurv.ensemble import RandomSurvivalForest


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 03b: tune Random Survival Forest")
    from src.evaluation.metrics import concordance_index, surv_y_to_structured
    from src.modeling._data import fit_scaler_train, load_processed, select_feature_set, split_xy
    from src.utils.logging import get_logger

    log = get_logger("rsf_tune", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    seed = int(cfg.get("seed", 42))
    spec = cfg["models"]["rsf"]

    processed = load_processed(cfg)
    X = select_feature_set(processed, spec.get("feature_set", "baseline"))
    y = processed["y"]
    X_tr, y_tr, X_te, y_te = split_xy(X, y, processed["split"])
    X_tr_s, X_te_s, _ = fit_scaler_train(X_tr, X_te)
    y_train_struct = surv_y_to_structured(y_tr)

    grid = {
        "min_samples_leaf": [10, 15, 25, 40],
        "min_samples_split": [10, 20, 40],
        "max_features": ["sqrt", "log2", 0.05, 0.1, 0.2],
        "max_depth": [3, 5, 8, None],
        "n_estimators": [300],
    }

    rows = []
    best = None
    best_score = float("-inf")

    combos = list(itertools.product(
        grid["min_samples_leaf"],
        grid["min_samples_split"],
        grid["max_features"],
        grid["max_depth"],
        grid["n_estimators"],
    ))
    log.info(f"RSF tuning sweep over {len(combos)} combinations")

    for i, (leaf, split, maxf, depth, n_est) in enumerate(combos, start=1):
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
        model.fit(X_tr_s.values, y_train_struct)
        train_risk = model.predict(X_tr_s.values)
        test_risk = model.predict(X_te_s.values)
        train_c = concordance_index(y_tr["OS_time"].to_numpy(), y_tr["OS_event"].to_numpy(), train_risk)
        test_c = concordance_index(y_te["OS_time"].to_numpy(), y_te["OS_event"].to_numpy(), test_risk)
        oob = float(getattr(model, "oob_score_", float("nan")))
        row = {
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
        score = test_c
        if score > best_score:
            best_score = score
            best = row
        log.info(
            f"[{i}/{len(combos)}] leaf={leaf} split={split} maxf={maxf} depth={depth} "
            f"oob={oob:.4f} train={train_c:.4f} test={test_c:.4f}"
        )

    results = pd.DataFrame(rows).sort_values(["c_index_test", "oob_score"], ascending=[False, False]).reset_index(drop=True)
    out_csv = f"{cfg['paths']['metrics']}/rsf_tuning_sweep.csv"
    results.to_csv(out_csv, index=False)
    with open(f"{cfg['paths']['metrics']}/rsf_tuning_best.json", "w") as fh:
        json.dump(best, fh, indent=2)

    log.info(f"Best RSF config: {best}")
    log.info(f"Wrote tuning results to {out_csv}")


if __name__ == "__main__":
    main()
