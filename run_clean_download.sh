#!/usr/bin/env bash
set -euo pipefail
RAW=/home/exe0x/mnt/datapool/tcga_survival/data/raw
find "$RAW" -maxdepth 1 \( \
  -name 'TCGA-LUAD' -o \
  -name 'TCGA-LUAD.clinical.tsv.gz' -o \
  -name 'TCGA-LUAD.survival.tsv.gz' -o \
  -name 'TCGA-LUAD.htseq_counts.tsv.gz' -o \
  -name 'TCGA-LUAD.gistic.tsv.gz' -o \
  -name 'TCGA-LUAD.methylation450.tsv.gz' -o \
  -name 'TCGA-LUAD.mirna.tsv.gz' -o \
  -name 'TCGA-LUAD.mutect2_snv.tsv.gz' -o \
  -name 'Wed_Apr_29_14_00_12_2026_*.tar.gz' -o \
  -name 'df.rds' -o -name 'results.rds' -o -name 'MANIFEST.txt' \
\) -print0 | xargs -0 rm -rf --
cd /home/exe0x/tcga_survival
exec /home/exe0x/tcga_survival/.venv/bin/python /home/exe0x/tcga_survival/tmp_resume_meth_screen.py >> /home/exe0x/mnt/datapool/tcga_survival/logs/download_gdc_clean_rerun.out 2>&1
