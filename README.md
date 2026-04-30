# TCGA Multi-Modal Survival Prediction

End-to-end pipeline for predicting **patient overall survival** from
TCGA cancer data using clinical, RNA-seq, somatic-mutation, copy-number,
DNA-methylation, and miRNA modalities, plus a modular **R Shiny atlas**
for exploring cross-cohort benchmark results, tumor-vs-normal biology,
unsupervised molecular structure, and survival patterns.

The project now has two linked outputs:

1. a **Python benchmark pipeline** that downloads, preprocesses, trains,
   and evaluates survival models per TCGA cohort; and
2. a **Shiny exploration layer** in `shiny_app/` that converts those
   saved outputs into an interactive computational oncology environment.

The default cohort is **TCGA-LUAD** (lung adenocarcinoma), but the
pipeline is structured so other TCGA cohorts can be swapped in by editing
`config.yaml` or by using the multi-cohort driver.

## What this repository contains

- Cox elastic-net, Random Survival Forest, DeepSurv, and MultiBranch
  survival modeling
- per-cohort train/test benchmark outputs and saved figures
- cross-cohort aggregation utilities
- a Shiny app for interactive exploration of benchmark and biology results
- documentation of exclusions and caveats, including the current
  `TCGA-LAML` ingestion failure

## What this repository does not contain

- raw TCGA downloads
- large processed matrices
- full trained model artefacts for every cohort
- Shiny app cache files

Those artefacts are intentionally regenerated from scripts or stored
outside Git history.

---

## Project layout

```
tcga_survival/
├── README.md
├── requirements.txt
├── config.yaml                 # all paths, URLs, hyperparameters
├── data/                       # optional local paths; defaults are set in config.yaml
│   ├── raw/                    # downloaded Xena/GDC files
│   └── processed/              # parquet matrices + train/test split
├── results/                    # optional local paths; defaults are set in config.yaml
│   ├── models/                 # serialized fitted models
│   ├── figures/                # KM curves, training plots, importances
│   └── metrics/                # per-model JSON + summary.csv
├── src/
│   ├── download/               # Xena (preferred) + GDC (stub)
│   ├── preprocessing/          # one module per modality + merge
│   ├── modeling/               # Cox / RSF / DeepSurv / MultiBranch
│   ├── evaluation/             # metrics, KM curves, importance
│   └── utils/                  # io, logging, seeding
└── scripts/                    # CLI entry points 01-06
```

---

## Modeling progression

| Step | Model | Lib | Inputs | Loss / Criterion |
|-----:|-------|-----|--------|------------------|
| 1 | Cox elastic-net | scikit-survival | concatenated features | partial log-likelihood, L1+L2 |
| 2 | Random Survival Forest | scikit-survival | concatenated features | log-rank splits |
| 3 | DeepSurv | PyTorch (custom loop) | concatenated features | Cox PH (Breslow) |
| 4 | MultiBranch | PyTorch (custom loop) | per-modality branches → fusion | Cox PH (Breslow) |

Every model is trained on the **same train/test split** and reports
**train + test C-index**, **median-risk Kaplan-Meier**, and a
**log-rank p-value**, so comparisons are honest.

---

## Setup

```bash
git clone <this-repo> && cd tcga_survival
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

The modeling pipeline only needs CPU PyTorch; GPU is auto-detected and
used if available.

For the Shiny app, create an R environment with the packages listed in
`shiny_app/README.md`.

---

## Common use cases

### 1. Train and evaluate one cohort

Use the six Python scripts in order with `config.yaml`.

### 2. Run a multi-cohort benchmark

Use `scripts/run_all_tcga_cohorts.py`, then aggregate with
`scripts/summarize_tcga_cohorts.py`.

### 3. Build and launch the Shiny app

Use the completed cohort outputs to build `shiny_app/data_cache/` via
`shiny_app/scripts/prepare_app_data.R`, validate with
`shiny_app/scripts/sanity_check.R`, then launch with
`shiny_app/scripts/run_app.R`.

See **Full rebuild order (benchmark + app)** below for the canonical
sequence.

---

## Running the pipeline

Each script is independent and reads `config.yaml` for paths/URLs/hyperparameters.
By default, raw data, processed matrices, logs, and model outputs are written under
`~/mnt/datapool/tcga_survival/` to avoid filling the home directory; edit the
`paths:` block in `config.yaml` to relocate them.

```bash
python scripts/01_download_data.py --config config.yaml
python scripts/02_preprocess_data.py --config config.yaml
python scripts/03_train_baselines.py --config config.yaml      # Cox + RSF
python scripts/04_train_deepsurv.py --config config.yaml
python scripts/05_train_multibranch_model.py --config config.yaml
python scripts/06_evaluate_models.py --config config.yaml      # summary table
```

Re-running step 01 is safe — it skips files already present in the configured
raw-data directory (`paths.data_raw`). Step 02 can be re-run after editing the
`preprocessing.*` block of `config.yaml`.

After step 06 you get:

* `results/metrics/summary.csv` — one row per model with C-index, KM
  log-rank p, and feature-set info.
* `results/figures/cindex_comparison.png` — side-by-side train vs test
  C-index bar chart.
* per-model KM plots, training curves, and feature importance
  visualizations next to it.

---

## Data source

UCSC Xena is used by default because it ships **harmonized OS and
days-to-event** for every TCGA cohort, which dramatically simplifies
clinical preprocessing. The full URL list is in `config.yaml ->
xena.urls`. The GDC API path (`src/download/download_gdc.py`) is
scaffolded as an opt-in fallback.

For TCGA-LUAD via the GDC Xena hub the seven files used are:

* clinical (`TCGA-LUAD.clinical.tsv.gz`)
* survival (`TCGA-LUAD.survival.tsv.gz`)
* RNA-seq FPKM (`TCGA-LUAD.htseq_fpkm.tsv.gz`)
* somatic mutations from MuTect2 (`TCGA-LUAD.mutect2_snv.tsv.gz`)
* GISTIC CNV (`TCGA-LUAD.gistic.tsv.gz`)
* Methylation 450K (`TCGA-LUAD.methylation450.tsv.gz`)
* miRNA expression (`TCGA-LUAD.mirna.tsv.gz`)

URLs occasionally change as Xena rebuilds the hub; if a download fails,
visit <https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Adenocarcinoma>
and update the URL in `config.yaml`.

---

## Switching cancer types

To run on a different TCGA cohort (BRCA, LIHC, GBM, …):

1. Set `project: TCGA-BRCA` in `config.yaml`.
2. Replace the seven URLs under `xena.urls.<modality>.url` with the
   matching cohort URLs from the GDC Xena hub.
3. Update each `filename:` to match.
4. Re-run from step 01.

No code change is required.

### Known cohort-specific GDC issue

In the multi-cohort GDC sweep, `TCGA-LAML` currently remains excluded because
`TCGAbiolinks` can fail on cohort-specific clinical metadata handling
(specifically around missing `disease_response` fields during `GDCquery_clinic`
and downstream `GDCprepare` calls). This is a downloader-layer incompatibility,
not a modeling result. The current aggregate summaries therefore omit `TCGA-LAML`
until that GDC/TCGAbiolinks path is hardened further.

---

## Adding a new modality

1. Add a preprocessing module under `src/preprocessing/` that exposes
   `def preprocess_<modality>(path, cfg) -> pd.DataFrame` returning a
   patient-indexed feature frame.
2. Wire its file URL into `config.yaml -> xena.urls.<modality>` and call
   it from `scripts/02_preprocess_data.py`.
3. Pass the frame into `merge_modalities(...)`. It will be saved as
   `data/processed/X_<modality>.parquet` and become available to the
   multi-branch model under `cfg.models.multibranch.modalities`.

The branch is built automatically as long as the parquet file exists;
if you also want it in the deep-baseline, append the prefix to
`X_multimodal` inside `merge_modalities.py`.

---

## Reproducibility

* `seed: 42` controls Python `random`, NumPy, and PyTorch RNGs (via
  `src/utils/seed.py`). Override per run by editing the config.
* All preprocessing decisions (gene counts, log transforms, CNV mode,
  variance thresholds) are driven by `config.yaml`.
* Every script mirrors stdout into `logs/<name>.log` so you have a
  trail of exactly what was done.
* All scaling, top-K feature selection, and imputation that depend on
  values are fit on **train only** inside the modeling scripts. Top-K
  variance pre-filtering at the preprocessing step is computed on the
  union but is unsupervised, so it does not leak survival information.
* Pre-trained model artifacts live under `results/models/`. For
  PyTorch models we save both the `state_dict` and a `*_meta.pkl`
  (scaler, feature names, in_features) so inference is reproducible.

---

## Full rebuild order (benchmark + app)

This repository now has two layers:

1. the **Python survival-modeling pipeline**, which creates per-cohort
   processed matrices, trained models, metrics, and figures; and
2. the **R Shiny app**, which reads those saved outputs and builds a
   lightweight `shiny_app/data_cache/` for interactive exploration.

If you are publishing this repository without raw data or large result
artefacts, the most important thing to provide is the exact rebuild order.
Use the steps below as the authoritative sequence.

### A. Quick path: run the app from existing cohort outputs

Use this if you already have a populated cohort-output tree such as
`~/mnt/datapool/tcga_survival_cohorts/`.

1. Clone the repository and create the Python + R environments.
2. Ensure the per-cohort outputs exist under a single root directory,
   with one subdirectory per cohort containing:
   - `data/raw/`
   - `data/processed/`
   - `results/metrics/`
   - `results/figures/`
3. Build the app cache:
   ```bash
   RSCRIPT=/path/to/Rscript
   export TCGA_SURVIVAL_COHORT_ROOT=/path/to/tcga_survival_cohorts
   $RSCRIPT shiny_app/scripts/prepare_app_data.R --no-methylation
   ```
4. Validate the cache:
   ```bash
   $RSCRIPT shiny_app/scripts/sanity_check.R
   ```
5. Launch the app:
   ```bash
   $RSCRIPT shiny_app/scripts/run_app.R
   ```

### B. Full path: rebuild the benchmark from scratch

#### Step 1 — environment setup

Create the Python environment for the modeling pipeline:

```bash
git clone <this-repo> && cd tcga_survival
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Create an R environment that includes the packages listed in
`shiny_app/README.md`. The app and cache-prep scripts require `Rscript`,
`data.table`, `jsonlite`, `stringr`, `shiny`, `bslib`, `ggplot2`, and
other listed packages; `arrow`/`pyarrow` support is needed at prep time.

#### Step 2 — choose output locations

For a single cohort, `config.yaml` controls all paths. For a multi-cohort
run, `scripts/run_all_tcga_cohorts.py` writes cohort-specific configs and
stores each cohort under a dedicated directory.

Expected multi-cohort layout:

```text
/path/to/tcga_survival_cohorts/
  TCGA-LUAD/
    config.yaml
    data/raw/
    data/processed/
    results/models/
    results/figures/
    results/metrics/
    logs/
    run_summary.json
  TCGA-BRCA/
  ...
  _summary/
```

#### Step 3 — build one cohort end-to-end

Before launching a full multi-cohort sweep, test the pipeline on one cohort
(first-choice: `TCGA-LUAD`) using the canonical script order:

```bash
python scripts/01_download_data.py --config config.yaml
python scripts/02_preprocess_data.py --config config.yaml
python scripts/03_train_baselines.py --config config.yaml
python scripts/04_train_deepsurv.py --config config.yaml
python scripts/05_train_multibranch_model.py --config config.yaml
python scripts/06_evaluate_models.py --config config.yaml
```

Outputs expected after step 06:
- `results/metrics/summary.csv`
- per-model metrics JSON files in `results/metrics/`
- saved figures in `results/figures/`
- trained model artefacts in `results/models/`

#### Step 4 — scale to many cohorts

After one cohort works, run the multi-cohort driver:

```bash
python scripts/run_all_tcga_cohorts.py
```

This script:
- iterates over the predefined TCGA cohort list,
- writes one cohort-specific `config.yaml` per cohort,
- runs steps 01–06 in order for each cohort,
- writes `run_summary.json` inside each cohort directory.

#### Step 5 — summarize across cohorts

After the cohort runs finish, aggregate the benchmark results:

```bash
python scripts/summarize_tcga_cohorts.py
```

Expected aggregate outputs:
- `.../_summary/aggregate_best_models.csv`
- `.../_summary/aggregate_best_models.json`
- `.../_summary/README.md`

#### Step 6 — build the Shiny app cache

Once the per-cohort results and aggregate summary exist, build the app cache:

```bash
RSCRIPT=/path/to/Rscript
export TCGA_SURVIVAL_COHORT_ROOT=/path/to/tcga_survival_cohorts
export TCGA_SURVIVAL_SUMMARY_DIR=/path/to/tcga_survival_cohorts/_summary
$RSCRIPT shiny_app/scripts/prepare_app_data.R --no-methylation
```

Optional flags:
- `--quick` to build only the lightest caches
- `--no-methylation` to skip very large methylation matrices
- `--cohorts=LUAD,BRCA,...` to restrict heavy cache generation
- `--no-umap` to skip precomputed UMAP

#### Step 7 — validate and launch the app

```bash
$RSCRIPT shiny_app/scripts/sanity_check.R
$RSCRIPT shiny_app/scripts/run_app.R
```

If launching on a remote server for laptop/browser access:

```bash
export HOST=0.0.0.0
export PORT=3838
$RSCRIPT shiny_app/scripts/run_app.R
```

Then open `http://<server-ip>:3838` from another machine on the same network.

## Rebuild table

| Step | Script | Main input | Main output |
|---|---|---|---|
| 1 | `scripts/01_download_data.py` | Xena/GDC URLs + config | raw cohort files |
| 2 | `scripts/02_preprocess_data.py` | raw cohort files | processed parquet matrices |
| 3 | `scripts/03_train_baselines.py` | processed matrices | Cox + RSF models/metrics |
| 4 | `scripts/04_train_deepsurv.py` | processed matrices | DeepSurv models/metrics |
| 5 | `scripts/05_train_multibranch_model.py` | processed matrices | MultiBranch models/metrics |
| 6 | `scripts/06_evaluate_models.py` | per-model metrics JSON | `summary.csv` + figures |
| 7 | `scripts/run_all_tcga_cohorts.py` | base config | full cohort tree |
| 8 | `scripts/summarize_tcga_cohorts.py` | cohort metrics | `_summary/aggregate_best_models.*` |
| 9 | `shiny_app/scripts/prepare_app_data.R` | cohort outputs + summary | `shiny_app/data_cache/` |
| 10 | `shiny_app/scripts/sanity_check.R` | app cache | validation report |
| 11 | `shiny_app/scripts/run_app.R` | app cache | running Shiny app |

## Important caveats for rebuilders

- `TCGA-LAML` is currently **excluded** due to a GDC/TCGAbiolinks clinical
  metadata incompatibility around missing fields such as `disease_response`.
  Treat this as a downloader/ingestion issue, not a benchmark result.
- Not every discovered TCGA cohort is complete. The app shows excluded /
  partial / incomplete status explicitly.
- Methylation is the heaviest modality by RAM and storage; skipping it is
  acceptable for a faster first rebuild.
- The Shiny app does **not** need the full raw-data tree at runtime; it
  primarily needs `shiny_app/data_cache/`, and optionally the cohort root if
  you want the saved Survival Modeling figures.
- For GitHub, prefer uploading the **scripts and documentation** needed to
  regenerate models rather than committing large trained model binaries.
  If exact frozen artefacts are important, publish them as release assets or
  in an external archive rather than in the main Git history.

---

## Implementation philosophy

* **Working > clever.** Cox/RSF baselines are the bar to clear, not the
  fallback after a deep model fails.
* **One split, one risk score, one evaluator.** Every model emits a
  patient-level risk and is graded with the same C-index + KM.
* **Modular branches.** Modalities can be missing without breaking the
  multi-branch model; the active branches are reported in the metrics.
* **Don't crown the deep model.** The summary table reports honest
  test-set C-index; if Cox elastic-net wins on this cohort, that's the
  result.
