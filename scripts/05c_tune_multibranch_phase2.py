from __future__ import annotations

import itertools
import json

import _bootstrap  # noqa: F401


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 05c: tune MultiBranch phase 2")
    from src.modeling import train_multibranch
    from src.utils.logging import get_logger

    log = get_logger("multibranch_phase2", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    base_spec = dict(cfg["models"]["multibranch"])
    base_spec["modalities"] = ["clinical", "rna", "mutation", "cnv", "methylation"]

    prefilters = [
        {"rna": 1000, "cnv": 1000, "methylation": 2000},
        {"rna": 500, "cnv": 500, "methylation": 1000},
        {"rna": 1500, "cnv": 1000, "methylation": 3000},
    ]
    branch_hidden_opts = [[64, 32], [32]]
    branch_embedding_opts = [16, 8]
    fusion_hidden_opts = [[64, 32], [32]]
    dropout_opts = [0.3, 0.5]
    weight_decay_opts = [1e-4, 1e-3]

    combos = list(itertools.product(
        prefilters,
        branch_hidden_opts,
        branch_embedding_opts,
        fusion_hidden_opts,
        dropout_opts,
        weight_decay_opts,
    ))
    log.info(f"MultiBranch phase2 tuning over {len(combos)} combinations")

    rows = []
    best = None
    best_score = float("-inf")

    for i, (pref, bh, be, fh, drop, wd) in enumerate(combos, start=1):
        cfg["models"]["multibranch"] = dict(base_spec)
        cfg["models"]["multibranch"]["branch_prefilter_top_k"] = pref
        cfg["models"]["multibranch"]["branch_hidden"] = bh
        cfg["models"]["multibranch"]["branch_embedding"] = be
        cfg["models"]["multibranch"]["fusion_hidden"] = fh
        cfg["models"]["multibranch"]["dropout"] = drop
        cfg["models"]["multibranch"]["weight_decay"] = wd
        cfg["models"]["multibranch"]["batch_norm"] = True
        log.info(f"[{i}/{len(combos)}] pref={pref} bh={bh} be={be} fh={fh} drop={drop} wd={wd}")
        out = train_multibranch(cfg)
        m = out["metrics"]
        row = {
            "branch_prefilter_top_k": pref,
            "branch_hidden": bh,
            "branch_embedding": be,
            "fusion_hidden": fh,
            "dropout": drop,
            "weight_decay": wd,
            "c_index_train": m["c_index_train"],
            "c_index_test": m["c_index_test"],
            "gap_train_test": m["c_index_train"] - m["c_index_test"],
            "best_val_cindex": m["best_val_cindex"],
            "epochs_trained": m["epochs_trained"],
        }
        rows.append(row)
        if row["c_index_test"] > best_score:
            best_score = row["c_index_test"]
            best = row
        log.info(
            f"MultiBranch phase2 result test={row['c_index_test']:.4f} train={row['c_index_train']:.4f} gap={row['gap_train_test']:.4f}"
        )

    with open(f"{cfg['paths']['metrics']}/multibranch_phase2_results.json", "w") as fh:
        json.dump(rows, fh, indent=2)
    with open(f"{cfg['paths']['metrics']}/multibranch_phase2_best.json", "w") as fh:
        json.dump(best, fh, indent=2)

    log.info(f"Best MultiBranch phase2 config: {best}")


if __name__ == "__main__":
    main()
