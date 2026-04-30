#!/usr/bin/env bash
set -euo pipefail
cd /home/exe0x/tcga_survival
exec /home/exe0x/tcga_survival/.venv/bin/python /home/exe0x/tcga_survival/tmp_resume_meth_screen.py >> /home/exe0x/mnt/datapool/tcga_survival/logs/download_gdc_resume_meth_screen.out 2>&1
