# TCGA Survival Atlas — R Shiny app

A polished, modular Shiny app that explores **TCGA pan-cancer survival benchmark
results** alongside **tumor-vs-normal molecular biology** with clustering and
survival interrogation. It is intended to feel like a guided computational
oncology environment, not a generic dashboard.

## Layout

```
shiny_app/
├── app.R                       Top-level Shiny entry point
├── global.R                    Pkg loads, paths, theme, eager cache reads
├── Dockerfile                  Container image (rocker/shiny base)
├── .dockerignore
├── R/
│   ├── utils_config.R          Path resolution + startup validation
│   ├── utils_io.R              File/cache loaders (lazy + safe)
│   ├── utils_data_prep.R       Reshaping, clustering, embeddings, markers
│   ├── utils_plots.R           Editorial ggplot/plotly helpers + plot export
│   ├── mod_overview.R          Page 1 — Overview
│   ├── mod_cohort_atlas.R      Page 2 — Cohort atlas
│   ├── mod_tumor_normal.R      Page 3 — Tumor vs normal biology
│   ├── mod_molecular_structure.R  Page 4 — Clustering & subgroups
│   ├── mod_survival_modeling.R Page 5 — Per-cohort survival
│   ├── mod_cross_cohort.R      Page 6 — Cross-cohort benchmark
│   └── mod_methods_caveats.R   Page 7 — Methods & caveats (incl. Diagnostics)
├── www/
│   ├── styles.css              Editorial styling
│   └── helpers.js              Small client-side niceties
├── scripts/
│   ├── prepare_app_data.R      Builds data_cache/ from existing pipeline outputs
│   ├── sanity_check.R          Headless validator (run before deploy / in CI)
│   ├── run_app.R               Tiny launcher honouring PORT/HOST env vars
│   └── _pq2csv.py              Tiny parquet→CSV helper used by the prep script
└── data_cache/                 Output of prepare_app_data.R (see MANIFEST.md)
```

## Prepare the data cache

The app loads only from `shiny_app/data_cache/`. Build it once with:

```bash
RSCRIPT=~/mnt/datapool/conda-envs/tcgabiolinks/bin/Rscript

# All cohorts, skip raw methylation (fast, ~3 min, ~38 MB cache)
$RSCRIPT shiny_app/scripts/prepare_app_data.R --no-methylation

# Add methylation for selected cohorts on demand
$RSCRIPT shiny_app/scripts/prepare_app_data.R --cohorts=LUAD
```

See `data_cache/MANIFEST.md` for what each cache file contains.

## Run the app

```bash
RSCRIPT=~/mnt/datapool/conda-envs/tcgabiolinks/bin/Rscript

# Quickest local launch
$RSCRIPT -e "shiny::runApp('~/tcga_survival/shiny_app', launch.browser = FALSE)"

# Or use the launcher (honours PORT / HOST / TCGA_SURVIVAL_* env vars)
$RSCRIPT shiny_app/scripts/run_app.R

# Headless sanity check (matches the in-app Diagnostics card; exits non-zero
# when any required cache is missing)
$RSCRIPT shiny_app/scripts/sanity_check.R
```

## Configuration via environment variables

All path assumptions are resolved through `R/utils_config.R::build_app_paths()`
which reads (in priority order) env vars → sensible defaults. No path is
hardcoded in any module.

| Variable                       | Default                                              | Purpose                                       |
| ------------------------------ | ---------------------------------------------------- | --------------------------------------------- |
| `TCGA_SURVIVAL_APP_ROOT`       | `~/tcga_survival/shiny_app`                          | App directory                                 |
| `TCGA_SURVIVAL_REPO_ROOT`      | parent of `APP_ROOT`                                 | Repo root (for README / METHODS rendering)    |
| `TCGA_SURVIVAL_CACHE_DIR`      | `$APP_ROOT/data_cache`                               | Output of `prepare_app_data.R`                |
| `TCGA_SURVIVAL_COHORT_ROOT`    | `~/mnt/datapool/tcga_survival_cohorts`               | Per-cohort raw + processed pipeline outputs   |
| `TCGA_SURVIVAL_SUMMARY_DIR`    | `$COHORT_ROOT/_summary`                              | Aggregate benchmark outputs                   |
| `TCGA_SURVIVAL_LOG_LEVEL`      | `info`                                               | `debug` / `info` / `warn` / `error`           |
| `TCGA_SURVIVAL_PYTHON`         | `~/venv/bin/python`                                  | Python with `pyarrow` for the parquet helper  |
| `PORT` / `HOST`                | `3838` / `0.0.0.0`                                   | Used by `scripts/run_app.R`                   |

## Deployment

### Local (developer machine)

```bash
$RSCRIPT shiny_app/scripts/sanity_check.R          # validate caches
$RSCRIPT -e "shiny::runApp('~/tcga_survival/shiny_app')"
```

### Posit Connect / Shiny Server (single-user)

1. Run `prepare_app_data.R` once locally so `data_cache/` is fully populated.
2. Push the entire `shiny_app/` directory **including `data_cache/`** to the
   server (the cache is roughly 122 MB at 18 cohorts).
3. Set environment variables in the deployment manifest so the app can find
   the cache and (optionally) the cohort root:
   ```
   TCGA_SURVIVAL_CACHE_DIR=/path/to/shiny_app/data_cache
   TCGA_SURVIVAL_COHORT_ROOT=/srv/research/tcga_survival_cohorts   # if mounted
   ```
4. The app **does not require** `TCGA_SURVIVAL_COHORT_ROOT` to be present —
   the Survival Modeling page will fall back to a calm "no saved figures"
   empty state, and the Methods → Diagnostics card will show the missing
   path explicitly. Everything else works from `data_cache/` alone.
5. The Python parquet helper is **only** needed at prep time. Deployment
   targets do not need Python or pyarrow.

### Containerized

A minimal `Dockerfile` is provided (rocker/shiny base):

```bash
cd shiny_app
docker build -t tcga-survival-atlas .

docker run --rm -p 3838:3838 \
    -v "$PWD/data_cache:/srv/cache:ro" \
    -e TCGA_SURVIVAL_CACHE_DIR=/srv/cache \
    tcga-survival-atlas
```

Mount the cohort root too (`-v /path:/srv/cohorts:ro -e TCGA_SURVIVAL_COHORT_ROOT=/srv/cohorts`)
if you want the saved figures on the Survival Modeling page. Otherwise that
page degrades gracefully.

The image's `HEALTHCHECK` runs `sanity_check.R` so orchestrators see an
unhealthy container when caches are misconfigured.

### Validation in CI

```yaml
- name: Validate Shiny cache
  run: Rscript shiny_app/scripts/sanity_check.R
```

Exit codes: `0` (all ok or warn-only), `1` (at least one missing artefact),
`2` (unhandled R error).

## Pages and what they answer

| Page                   | Scientific question                                                  |
| ---------------------- | -------------------------------------------------------------------- |
| Overview               | What is this project, what are the headline findings?                |
| Cohort atlas           | What data exist for each cohort, who's complete vs excluded?         |
| Tumor vs normal        | How do molecular values differ between tumor and normal?             |
| Molecular structure    | What sample subtypes emerge — and how do they relate to clinical/survival? |
| Survival modeling      | How well do models predict survival within a cohort?                 |
| Cross-cohort benchmark | What patterns hold (or don't) across cohorts?                        |
| Methods                | Methodology, caveats, exclusions (esp. LAML).                        |

## Package requirements

R packages are listed in `global.R`. The full set:

* Core: `shiny`, `bslib`, `htmltools`, `bsicons`, `thematic`, `shinyWidgets`,
  `DT`, `reactable`
* Plotting: `ggplot2`, `plotly`, `pheatmap`, `patchwork`, `cowplot`, `viridis`,
  `RColorBrewer`, `ggrepel`, `scales`, `base64enc`
* Data: `data.table`, `dplyr`, `tidyr`, `stringr`, `purrr`, `readr`,
  `jsonlite`, `arrow`, `yaml`, `glue`, `fs`
* Stats / clustering: `cluster`, `factoextra`, `uwot`, `matrixStats`
* Survival: `survival`, `survminer`

The app degrades gracefully when optional packages are missing: e.g. if
`plotly` isn't available, plots are rendered as static `ggplot` instead of
interactive widgets; if `DT` is missing, plain `tableOutput` is used.

## Defaults

* Cohort: `TCGA-LUAD`
* Modality: `RNA`
* Clustering: hierarchical (Ward.D2), k = 3, scaled features
* Embedding: PCA (uses precomputed scores from cache when available)

## Known limitations

* `TCGA-LAML` is excluded — TCGAbiolinks could not pull required clinical
  fields (notably `disease_response`); shown explicitly in the app rather
  than dropped silently.
* `TCGA-OV` (and any other not-yet-finished cohort) is shown as
  `incomplete`. Light caches still include it; heavy caches skip it.
* CNV (GISTIC) and somatic mutation (MAF) only contain tumor samples, so
  tumor-vs-normal feature summaries / volcanos are intentionally not
  produced for those modalities. Some normal-dependent views may show a
  calm "no normals available" empty state for these modalities.
* Methylation raw matrices are large (~1–3 GB / cohort); the default
  prep run includes them but `--no-methylation` is available for fast
  rebuilds.
* C-index measures discrimination only — not calibration.
* Test sets are small for several cohorts; treat differences below
  ~0.05 in C-index as inside the noise floor.
* The Survival Modeling page renders **precomputed** model figures from
  the cohort root — without that mount the page degrades to "no saved
  figures" but everything else still works from `data_cache/`.

### Remaining blockers for a hardened production deployment

These are intentionally out of scope for this iteration but worth noting
for any next deployment pass:

* No automated regression tests — `sanity_check.R` only validates
  presence/shape of caches, not numerical equivalence between runs.
* No authentication layer; the app assumes the host (Posit Connect /
  reverse proxy) handles auth.
* Bookmarkable URLs encode UI inputs, not server-side derived state, so
  shared links are stable across cache refreshes but **not** across
  schema changes (e.g., re-running prep with different feature counts).
* The clustering page recomputes on every parameter change. For very
  large cohorts × modalities this is acceptable on a developer machine
  but may need `bindCache()` / debouncing under heavy concurrent load.
* `shiny::reactivePoll()` is not used to watch the cache directory, so
  refreshing the underlying data requires restarting the app.

## UI / UX design

* **Three-pane layout** on every analytical page: left rail = controls,
  center = primary evidence/plots, right rail = interpretation, caveats
  and a "current selection" summary. The right rail is not decorative —
  it is the page's reading guide.
* **Editorial page headers** — every page opens with an eyebrow line
  ("modality view · per cohort"), a strong title and a 1–2 sentence
  lede that frames the scientific question.
* **Centralized visual language** — all categorical colours come from
  `TCGA_PALETTE` in `global.R` and the matched CSS variables in
  `www/styles.css`. Tumor vs normal, cluster, model family and stage
  use the same colours wherever they appear.
* **Dynamic interpretation** — `interp_tumor_normal()`,
  `interp_clustering()` and `interp_survival_cohort()` generate
  context-aware sentences (e.g. *"Tumor–normal separation is visually
  strong: 312 of 1000 features cross the |Δ|≥1 / FDR<1% threshold"*,
  *"Survival differs strongly between clusters (logrank p = 0.003)"*,
  *"This cohort has 22 test events — interpret survival separations
  cautiously; CIs will be wide"*).
* **Calm empty states** — every plot/table has a structured
  `empty_state(title, body, icon)` panel for the *no normals*, *cohort
  excluded*, *missing cache* cases instead of raw red errors.
* **Status pills** — `status_pill()` renders `completed` / `partial` /
  `incomplete` / `excluded` as Bootstrap-coloured badges; the legend
  lives on the Cohort Atlas page.
* **Where to go next** — the Overview page ends with four clickable
  cards that take the user to the Cohort Atlas, Tumor-vs-Normal,
  Molecular Structure or Cross-Cohort views (handled by a small
  `Shiny.setInputValue('._goto', …)` ↔ `nav_select()` bridge).
* **Loading indicators** — clustering and on-the-fly UMAP show a
  `withProgress()` overlay so the user knows the page is computing
  rather than frozen.
* **Kaplan–Meier panels** include a log-rank p annotation in the
  bottom-left, making the figure self-contained for screenshots.
* **Discoverable but secondary downloads** — every page has a single
  pinned `Download …` button in the left rail with a consistent label.
  CSS keeps these buttons quiet (outline-secondary, not primary).
* **No business-dashboard chrome** — minty + Inter + Source Sans 3,
  light background, restrained accents, generous typography, no large
  KPI tiles, subtle card borders.

## Features

### Tumor vs normal
* Box / violin / jitter for any selected feature
* Sortable, searchable top-differential table with Welch's t-test
  `pvalue` and BH-corrected `fdr`
* Volcano plot (`diff_mean` vs `−log10 p`) with auto-tuned thresholds
  and `ggrepel` labelling of top hits
* Top-features × samples heatmap (ComplexHeatmap; pheatmap / base
  fallback) annotated by `sample_type`, `stage`, `sex`
* PCA scatter colored tumor vs normal (uses precomputed PCA cache)
* Mutation frequency view appears when modality = mutation

### Molecular structure
* Hierarchical (Ward.D2) or k-means clustering with adjustable k and
  feature-scaling toggle
* PCA *or* UMAP embedding — both precomputed in the data-prep step,
  loaded instantly from cache
* Color points by cluster / sample_type / stage / sex / OS_event
* Cluster-size bars + cluster vs clinical stacked bar
* Kaplan–Meier curves by cluster (when survival labels are available)
* Top features per cluster (one-vs-rest mean difference) shown as
  table + cluster × marker heatmap

### Reproducibility / sharing
* URL bookmarking — every input (cohort, modality, k, embedding,
  filters, …) is encoded in the URL via `enableBookmarking = "url"`.
  Click the **bookmark** link in the top navbar to grab a shareable
  URL of the current view.
* "Current solution" / "At a glance" right-rail panels show the active
  cohort/modality/scope/method/k/embedding/color-by, so the analysis
  is self-describing in screenshots.
* CSV downloads for benchmark, cohort metadata, tumor/normal
  rankings, cohort summary, and cluster assignments.
* PNG plot exports for the volcano (Tumor vs Normal) and Kaplan–Meier
  (Molecular Structure) panels via `download_plot_handler()`.

### Validation & runtime safety
* `R/utils_config.R::validate_app_data()` runs at startup and powers
  the Methods → **Diagnostics** card. Status pills surface
  ok / warn / missing items so a misconfigured cache is obvious to
  the user, not buried as a stack trace.
* `scripts/sanity_check.R` runs the same checks headlessly — useful
  for CI gates and the Docker `HEALTHCHECK`.
* Missing optional artefacts (figures, normals, modality caches, even
  the cohort root) degrade to calm `empty_state(...)` panels rather
  than red Shiny errors.

## Reproducibility / sessions

Bookmarking is enabled (`enableBookmarking = "url"`) — the navbar has a
**bookmark** link that pops up a copyable, self-contained URL encoding
the entire UI state. Internal table-state and download triggers are
excluded from the URL via `setBookmarkExclude()` to keep links short
and stable.
