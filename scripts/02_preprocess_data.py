"""Run every modality preprocessor and write data/processed/*.parquet."""
from __future__ import annotations

from pathlib import Path

import _bootstrap  # noqa: F401


def _resolve(cfg, name: str):
    """Resolve a modality file in data/raw for either Xena or GDC naming."""
    raw_dir = Path(cfg["paths"]["data_raw"])

    # Prefer configured Xena filename when present.
    entries = cfg.get("xena", {}).get("urls", {})
    entry = entries.get(name)
    if entry:
        p = raw_dir / entry["filename"]
        if p.exists():
            return p

    project = cfg.get("gdc", {}).get("project_id") or cfg.get("project", "TCGA-LUAD")
    gdc_fallbacks = {
        "clinical": [f"{project}.clinical.tsv.gz"],
        "survival": [f"{project}.survival.tsv.gz"],
        "rnaseq": [f"{project}.htseq_counts.tsv.gz", f"{project}.htseq_fpkm.tsv.gz"],
        "mutation": [f"{project}.mutect2_snv.tsv.gz"],
        "cnv": [f"{project}.gistic.tsv.gz"],
        "methylation": [f"{project}.methylation450.tsv.gz"],
        "mirna": [f"{project}.mirna.tsv.gz"],
    }
    for fn in gdc_fallbacks.get(name, []):
        p = raw_dir / fn
        if p.exists():
            return p
    return None


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 02: preprocess TCGA modalities")
    from src.preprocessing import (
        merge_modalities,
        preprocess_clinical,
        preprocess_cnv,
        preprocess_methylation,
        preprocess_mirna,
        preprocess_mutations,
        preprocess_rnaseq,
    )
    from src.utils.logging import get_logger

    log = get_logger("script_preprocess", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))

    clinical_path = _resolve(cfg, "clinical")
    survival_path = _resolve(cfg, "survival")
    if clinical_path is None:
        raise FileNotFoundError("Clinical file missing in data/raw — run scripts/01_download_data.py first.")

    clinical = preprocess_clinical(clinical_path, survival_path, cfg)

    rna = mut = cnv = meth = mir = None
    if cfg.get("preprocessing", {}).get("rnaseq", {}).get("enabled", True):
        if (p := _resolve(cfg, "rnaseq")) is not None:
            rna = preprocess_rnaseq(p, cfg)
        else:
            log.warning("RNA-seq file not found — continuing without it.")
    else:
        log.info("RNA-seq preprocessing disabled in config; skipping.")

    if cfg.get("preprocessing", {}).get("mutations", {}).get("enabled", True):
        if (p := _resolve(cfg, "mutation")) is not None:
            mut = preprocess_mutations(p, cfg)
        else:
            log.warning("Mutation file not found — continuing without it.")
    else:
        log.info("Mutation preprocessing disabled in config; skipping.")

    if cfg.get("preprocessing", {}).get("cnv", {}).get("enabled", True):
        if (p := _resolve(cfg, "cnv")) is not None:
            cnv = preprocess_cnv(p, cfg)
        else:
            log.warning("CNV file not found — continuing without it.")
    else:
        log.info("CNV preprocessing disabled in config; skipping.")

    if cfg.get("preprocessing", {}).get("methylation", {}).get("enabled", True):
        if (p := _resolve(cfg, "methylation")) is not None:
            meth = preprocess_methylation(p, cfg)
        else:
            log.warning("Methylation file not found — continuing without it.")
    else:
        log.info("Methylation preprocessing disabled in config; skipping.")

    if cfg.get("preprocessing", {}).get("mirna", {}).get("enabled", True):
        if (p := _resolve(cfg, "mirna")) is not None:
            mir = preprocess_mirna(p, cfg)
        else:
            log.warning("miRNA file not found — continuing without it.")
    else:
        log.info("miRNA preprocessing disabled in config; skipping.")

    merge_modalities(clinical, rna, mut, cnv, meth, mir, cfg)
    log.info("Preprocessing complete.")


if __name__ == "__main__":
    main()
