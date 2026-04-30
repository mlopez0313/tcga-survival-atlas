from __future__ import annotations

import copy
import json
import subprocess
from pathlib import Path

import yaml

TCGA_PROJECTS = [
    "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL", "TCGA-COAD",
    "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC", "TCGA-KICH", "TCGA-KIRC",
    "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
    "TCGA-MESO", "TCGA-OV", "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ",
    "TCGA-SARC", "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
    "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM",
]

REPO = Path("/home/exe0x/tcga_survival")
BASE_CONFIG = REPO / "config.yaml"
OUT_ROOT = Path("/home/exe0x/mnt/datapool/tcga_survival_cohorts")
PYTHON = REPO / ".venv/bin/python"


def cohort_cfg(base: dict, project: str) -> dict:
    cfg = copy.deepcopy(base)
    slug = project
    cfg["project"] = project
    cfg.setdefault("gdc", {})["project_id"] = project
    cfg["paths"] = {
        "data_raw": str(OUT_ROOT / slug / "data/raw"),
        "data_processed": str(OUT_ROOT / slug / "data/processed"),
        "results": str(OUT_ROOT / slug / "results"),
        "models": str(OUT_ROOT / slug / "results/models"),
        "figures": str(OUT_ROOT / slug / "results/figures"),
        "metrics": str(OUT_ROOT / slug / "results/metrics"),
        "logs": str(OUT_ROOT / slug / "logs"),
    }
    if "xena" in cfg and "urls" in cfg["xena"]:
        for name, entry in cfg["xena"]["urls"].items():
            if isinstance(entry, dict):
                for key in ["url", "filename"]:
                    if key in entry and isinstance(entry[key], str):
                        entry[key] = entry[key].replace("TCGA-LUAD", project)
    return cfg


def run_step(config_path: Path, script_name: str, log_path: Path) -> int:
    cmd = [str(PYTHON), str(REPO / "scripts" / script_name), "--config", str(config_path)]
    with open(log_path, "w") as fh:
        proc = subprocess.run(cmd, cwd=REPO, stdout=fh, stderr=subprocess.STDOUT)
    return proc.returncode


def main() -> None:
    base = yaml.safe_load(BASE_CONFIG.read_text())
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary = []

    for project in TCGA_PROJECTS:
        cohort_dir = OUT_ROOT / project
        cohort_dir.mkdir(parents=True, exist_ok=True)
        cfg = cohort_cfg(base, project)
        cfg_path = cohort_dir / "config.yaml"
        cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))

        logs_dir = cohort_dir / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)

        steps = [
            ("01_download_data.py", logs_dir / "01_download_data.out"),
            ("02_preprocess_data.py", logs_dir / "02_preprocess_data.out"),
            ("03_train_baselines.py", logs_dir / "03_train_baselines.out"),
            ("04_train_deepsurv.py", logs_dir / "04_train_deepsurv.out"),
            ("05_train_multibranch_model.py", logs_dir / "05_train_multibranch_model.out"),
            ("06_evaluate_models.py", logs_dir / "06_evaluate_models.out"),
        ]

        status = {"project": project}
        for script_name, log_path in steps:
            rc = run_step(cfg_path, script_name, log_path)
            status[script_name] = rc
            if rc != 0:
                status["failed_at"] = script_name
                break
        summary.append(status)
        (cohort_dir / "run_summary.json").write_text(json.dumps(status, indent=2))

    (OUT_ROOT / "all_cohorts_summary.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
