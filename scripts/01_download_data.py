"""Download TCGA modality files defined in config.yaml -> data/raw/."""
from __future__ import annotations

import _bootstrap  # noqa: F401  (sys.path side effect)


def main() -> None:
    cfg, _ = _bootstrap.setup("Step 01: download TCGA data")
    source = cfg.get("data_source", "xena").lower()
    if source == "xena":
        from src.download.download_xena import download_xena_assets
        download_xena_assets(cfg)
    elif source == "gdc":
        from src.download.download_gdc import download_gdc_assets
        download_gdc_assets(cfg)
    else:
        raise ValueError(f"Unknown data_source: {source}")


if __name__ == "__main__":
    main()
