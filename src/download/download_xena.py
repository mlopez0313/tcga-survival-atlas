"""Download TCGA modality files from UCSC Xena.

URLs are read from `config.yaml -> xena.urls`, so swapping cohorts only
requires editing the YAML, not this code.

Files are streamed to `data/raw/`. Already-present non-empty files are
skipped, which makes the script idempotent and re-runnable.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import requests
from tqdm import tqdm

from ..utils.logging import get_logger


def _download_one(url: str, dest: Path, log) -> bool:
    """Stream-download a single file. Returns True on success / skip."""
    if dest.exists() and dest.stat().st_size > 0:
        log.info(f"[skip] {dest.name} already present ({dest.stat().st_size:,} bytes)")
        return True

    log.info(f"[get ] {url}  ->  {dest}")
    try:
        with requests.get(url, stream=True, timeout=60) as resp:
            resp.raise_for_status()
            total = int(resp.headers.get("Content-Length", 0))
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.parent.mkdir(parents=True, exist_ok=True)
            with open(tmp, "wb") as fh, tqdm(
                total=total or None,
                unit="B",
                unit_scale=True,
                desc=dest.name,
                leave=False,
            ) as pbar:
                for chunk in resp.iter_content(chunk_size=1 << 15):
                    if not chunk:
                        continue
                    fh.write(chunk)
                    pbar.update(len(chunk))
            tmp.replace(dest)
        return True
    except Exception as exc:  # noqa: BLE001 - we want to keep going
        log.error(f"failed: {url}  ({exc})")
        if dest.exists() and dest.stat().st_size == 0:
            dest.unlink(missing_ok=True)
        return False


def download_xena_assets(cfg: Dict[str, Any]) -> Dict[str, Path]:
    """Download every entry in cfg['xena']['urls']. Returns name -> path."""
    log = get_logger("download_xena", log_dir=cfg["paths"]["logs"], level=cfg.get("log_level", "INFO"))
    raw_dir = Path(cfg["paths"]["data_raw"])
    raw_dir.mkdir(parents=True, exist_ok=True)

    urls: Dict[str, Dict[str, str]] = cfg["xena"]["urls"]
    out: Dict[str, Path] = {}

    for name, entry in urls.items():
        url = entry["url"]
        filename = entry.get("filename") or url.split("/")[-1]
        dest = raw_dir / filename
        if _download_one(url, dest, log):
            out[name] = dest

    log.info(f"Downloaded / verified {len(out)}/{len(urls)} modalities into {raw_dir}")
    return out
