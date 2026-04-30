from __future__ import annotations

import json

import _bootstrap  # noqa: F401


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 05b: ablate MultiBranch modality sets")
    from src.modeling import train_multibranch
    from src.utils.logging import get_logger

    log = get_logger("multibranch_ablate", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    base_spec = dict(cfg["models"]["multibranch"])

    modality_sets = [
        ["clinical", "rna"],
        ["clinical", "rna", "mutation"],
        ["clinical", "rna", "cnv"],
        ["clinical", "rna", "mirna"],
        ["clinical", "rna", "mutation", "cnv", "mirna"],  # full minus methylation
        ["clinical", "rna", "mutation", "methylation", "mirna"],  # full minus cnv
        ["clinical", "rna", "mutation", "cnv", "methylation"],  # full minus mirna
        ["clinical", "rna", "mutation", "cnv", "methylation", "mirna"],
    ]

    rows = []
    best = None
    best_score = float("-inf")

    for mods in modality_sets:
        cfg["models"]["multibranch"] = dict(base_spec)
        cfg["models"]["multibranch"]["modalities"] = mods
        log.info(f"Running MultiBranch ablation with modalities={mods}")
        out = train_multibranch(cfg)
        m = out["metrics"]
        row = {
            "modalities": mods,
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
            f"MultiBranch mods={mods} train={row['c_index_train']:.4f} "
            f"test={row['c_index_test']:.4f} gap={row['gap_train_test']:.4f}"
        )

    with open(f"{cfg['paths']['metrics']}/multibranch_ablation_results.json", "w") as fh:
        json.dump(rows, fh, indent=2)
    with open(f"{cfg['paths']['metrics']}/multibranch_ablation_best.json", "w") as fh:
        json.dump(best, fh, indent=2)

    log.info(f"Best MultiBranch ablation: {best}")


if __name__ == "__main__":
    main()
