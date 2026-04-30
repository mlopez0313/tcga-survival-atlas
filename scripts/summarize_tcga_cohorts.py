from __future__ import annotations

import csv
import json
from pathlib import Path

OUT_ROOT = Path('/home/exe0x/mnt/datapool/tcga_survival_cohorts')
OUT_DIR = OUT_ROOT / '_summary'
OUT_DIR.mkdir(parents=True, exist_ok=True)

rows = []
for summary_csv in sorted(OUT_ROOT.glob('TCGA-*/results/metrics/summary.csv')):
    cohort = summary_csv.parts[-4]
    with open(summary_csv, newline='') as fh:
        data = list(csv.DictReader(fh))
    usable = []
    for r in data:
        try:
            c = float(r['c_index_test'])
        except Exception:
            continue
        usable.append((c, r))
    if not usable:
        continue
    usable.sort(key=lambda x: x[0], reverse=True)
    best = usable[0][1]
    rows.append({
        'cohort': cohort,
        'best_model': best.get('model', ''),
        'test_c_index': best.get('c_index_test', ''),
        'n_train': best.get('n_train', ''),
        'n_test': best.get('n_test', ''),
        'events_train': best.get('events_train', ''),
        'events_test': best.get('events_test', ''),
    })

rows.sort(key=lambda r: r['cohort'])

csv_path = OUT_DIR / 'aggregate_best_models.csv'
with open(csv_path, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else [
        'cohort', 'best_model', 'test_c_index', 'n_train', 'n_test', 'events_train', 'events_test'
    ])
    w.writeheader()
    w.writerows(rows)

json_path = OUT_DIR / 'aggregate_best_models.json'
json_path.write_text(json.dumps(rows, indent=2))

md_path = OUT_DIR / 'README.md'
with open(md_path, 'w') as fh:
    fh.write('# TCGA cohort aggregate summary\n\n')
    fh.write('| Cohort | Best model | Test C-index | n_train | n_test | events_train | events_test |\n')
    fh.write('|---|---|---:|---:|---:|---:|---:|\n')
    for r in rows:
        fh.write(f"| {r['cohort']} | {r['best_model']} | {r['test_c_index']} | {r['n_train']} | {r['n_test']} | {r['events_train']} | {r['events_test']} |\n")

print(csv_path)
print(json_path)
print(md_path)
