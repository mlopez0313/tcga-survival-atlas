import os
import sys
from pathlib import Path

ENV = Path('/home/exe0x/mnt/datapool/conda-envs/tcgabiolinks')
os.environ['PATH'] = f"{ENV / 'bin'}:{os.environ.get('PATH','')}"
os.environ['R_LIBS_USER'] = '/tmp/tcgabiolinks-rlib'

repo = Path('/home/exe0x/tcga_survival')
sys.path.insert(0, str(repo))
sys.path.insert(0, str(repo / 'scripts'))
import _bootstrap
from src.download.download_gdc import download_gdc_assets

cfg, _ = _bootstrap.setup('resume meth')
out = download_gdc_assets(cfg)
print({k: str(v) for k, v in out.items()})
