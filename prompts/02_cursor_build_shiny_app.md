Build a polished **R Shiny app** for the project in `~/tcga_survival` that explores both:

1. **TCGA multi-cohort survival benchmark results**
2. **Tumor vs normal molecular data across modalities** with clustering and survival interrogation

The app should **not** feel like a generic dashboard. It should feel like a **guided computational oncology exploration environment** organized around scientific questions.

# Project context

Relevant project roots:
- code repo: `~/tcga_survival`
- main outputs/data: `~/mnt/datapool/tcga_survival`
- pan-cohort outputs: `~/mnt/datapool/tcga_survival_cohorts`

Important existing files:
- aggregate cohort summary:
  - `~/mnt/datapool/tcga_survival_cohorts/_summary/aggregate_best_models.csv`
  - `~/mnt/datapool/tcga_survival_cohorts/_summary/aggregate_best_models.json`
  - `~/mnt/datapool/tcga_survival_cohorts/_summary/README.md`
- manuscript/results summary:
  - `~/tcga_survival/RESULTS_SUMMARY.md`
- completed cohort result dirs:
  - `~/mnt/datapool/tcga_survival_cohorts/TCGA-*/results/metrics/summary.csv`
  - `~/mnt/datapool/tcga_survival_cohorts/TCGA-*/results/figures/*`
  - `~/mnt/datapool/tcga_survival_cohorts/TCGA-*/results/models/*`
- per-cohort processed matrices:
  - `~/mnt/datapool/tcga_survival_cohorts/TCGA-*/data/processed/*.parquet`
- raw cohort files:
  - `~/mnt/datapool/tcga_survival_cohorts/TCGA-*/data/raw/*`
- app-specific caches, if prepared first:
  - `~/tcga_survival/shiny_app/data_cache/*`

Known data note:
- `TCGA-LAML` is currently excluded due to a GDC/TCGAbiolinks ingestion incompatibility involving missing clinical metadata fields like `disease_response`
- this should be explicitly surfaced in the app as an omitted cohort / caveat, not as missing silently

# App goals

The app must communicate:

1. **Cohort landscape**
   - what cohorts were analyzed
   - tumor/normal counts if available
   - event/censor counts
   - modality availability
   - completed vs excluded cohorts

2. **Survival benchmark results**
   - cross-cohort model comparison
   - cohort-specific model performance
   - KM curves and saved benchmark figures
   - no universal winner across cohorts
   - classical baselines remain competitive

3. **Tumor vs normal molecular biology**
   - distributions of RNA / miRNA / methylation / CNV values
   - per-feature tumor vs normal comparisons
   - multivariate structure like PCA/UMAP
   - clustering and subgroup discovery

4. **Clustering and subgroup interrogation**
   - clustering by clinical labels like stage
   - unsupervised clustering by modality
   - survival comparison by cluster
   - expression/methylation/etc. comparisons by cluster

# Design principles

- Do **not** build a generic KPI dashboard
- Organize around scientific questions, not file types
- Use a left control rail, central visualization area, right interpretation/caveat panel where appropriate
- Use spacious, publication-like layout
- Prefer high-value scientific plots over decorative widgets
- Avoid pie charts
- Add interpretation/help text directly in the app
- Use sensible defaults:
  - default cohort: `TCGA-LUAD`
  - default modality: RNA
  - default clustering: hierarchical / top variable features
- Make navigation intuitive and non-generic

# Required app structure

Create the app under:

- `~/tcga_survival/shiny_app/`

Use this structure:

- `shiny_app/app.R`
- `shiny_app/global.R`
- `shiny_app/R/`
- `shiny_app/R/mod_overview.R`
- `shiny_app/R/mod_cohort_atlas.R`
- `shiny_app/R/mod_tumor_normal.R`
- `shiny_app/R/mod_molecular_structure.R`
- `shiny_app/R/mod_survival_modeling.R`
- `shiny_app/R/mod_cross_cohort.R`
- `shiny_app/R/mod_methods_caveats.R`
- `shiny_app/R/utils_io.R`
- `shiny_app/R/utils_plots.R`
- `shiny_app/R/utils_data_prep.R`
- `shiny_app/www/`
- `shiny_app/www/styles.css`
- `shiny_app/www/helpers.js`
- `shiny_app/data_cache/`
- `shiny_app/README.md`

If needed, also create:
- `shiny_app/scripts/prepare_app_data.R`

# Preferred R packages

Use these unless there is a strong reason not to:

Core app:
- `shiny`
- `bslib`
- `htmltools`
- `htmlwidgets`
- `thematic`
- `shinyWidgets`
- `DT`
- `reactable`
- `crosstalk` optional

Plotting:
- `ggplot2`
- `plotly`
- `ComplexHeatmap` or `pheatmap`
- `patchwork`
- `cowplot`
- `viridis`
- `RColorBrewer`
- `ggrepel`

Data:
- `data.table`
- `dplyr`
- `tidyr`
- `stringr`
- `purrr`
- `readr`
- `jsonlite`
- `arrow`
- `yaml`
- `glue`
- `fs`

Dimensionality reduction / clustering:
- `stats`
- `cluster`
- `factoextra`
- `uwot`
- `Rtsne` (optional)
- `matrixStats`
- `proxy` (optional)
- `NMF` optional, only if easy
- `ConsensusClusterPlus` optional, not required in v1

Survival:
- `survival`
- `survminer`

If there are image display needs for saved PNGs:
- use standard HTML tags / `renderImage`

# Data/backend expectations

The app should be built around **precomputed results** where possible.

Prefer using app-specific caches if present in:
- `~/tcga_survival/shiny_app/data_cache/`

If caches are missing, gracefully fall back to direct reads from existing summary and processed files where feasible.

Use these sources:

## Benchmark summaries
- `_summary/aggregate_best_models.csv`
- per-cohort `results/metrics/summary.csv`

## Saved benchmark figures
- per-cohort `results/figures/*.png`

## Processed feature matrices
- per-cohort `data/processed/X_*.parquet`
- per-cohort `data/processed/y_survival.parquet`

## Raw molecular files
Use raw files only if needed to preserve tumor/normal sample structure and if caches are absent.

# App pages and exact intent

## 1. Overview
Scientific question:
- What is this project and what are the main findings?

Must include:
- title banner
- short description
- cards for:
  - completed cohorts
  - excluded cohorts
  - modalities
  - best-performing model family by win count
- cross-cohort performance heatmap
- win-count bar plot
- concise conclusion text
- explicit note that `TCGA-LAML` is omitted for ingestion reasons

## 2. Cohort Atlas
Scientific question:
- What data exist for each cohort?

Must include:
- cohort metadata table
- completion status
- modality availability heatmap
- sample/event summary plots
- if possible tumor vs normal counts by cohort
- filters:
  - completed only
  - min sample count
  - modality present
- right-side cohort details panel

## 3. Tumor–Normal Biology
Scientific question:
- How do molecular values differ between tumor and normal samples?

Must include:
- controls:
  - cohort
  - modality
  - feature selection/search
  - transform option
- plots:
  - boxplot / violin / jitter for selected feature
  - density/ridge optional
  - top feature heatmap
  - PCA / UMAP colored by tumor vs normal
- summary stats:
  - n tumor
  - n normal
  - mean/median by group
  - missingness
- support at least:
  - RNA
  - miRNA
  - methylation
  - CNV
- mutation can be a separate frequency-oriented view if helpful

## 4. Molecular Structure
Scientific question:
- What sample structure or subtypes emerge, and how do they relate to clinical variables and survival?

Must include:
- controls:
  - cohort
  - modality
  - clustering method (`hierarchical`, `kmeans`; optional more)
  - number of clusters
  - feature subset mode (`top variable`, maybe user-selected)
  - dimensionality reduction method (`PCA`, `UMAP`)
  - color-by variable (`cluster`, `stage`, `tumor_normal`, `sex`, `event`)
- plots:
  - embedding scatter
  - clustered heatmap
  - cluster size bar plot
  - stacked bar plot of cluster vs clinical variable
  - Kaplan–Meier curves by cluster
- marker table:
  - top features differing across clusters
- right-side interpretation panel summarizing selected cluster solution

Important:
- support grouping/interrogation by sensible clinical variables like stage
- support unsupervised clustering based on RNA, miRNA, methylation, CNV
- allow downstream interrogation of survival and expression based on cluster labels

## 5. Survival Modeling
Scientific question:
- How well do models predict survival within a cohort?

Must include:
- controls:
  - cohort
  - model(s)
- outputs:
  - per-cohort summary table from `summary.csv`
  - train vs test metric comparison
  - C-index comparison plot
  - saved KM curves
  - saved training curves and feature importance/coefficients if available
- for selected cohort, show all available saved figures from `results/figures/`
- emphasize caveats for small cohorts / low event counts

## 6. Cross-Cohort Benchmark
Scientific question:
- What patterns hold across cohorts?

Must include:
- heatmap of cohort × model test C-index
- win count bar plot
- scatter plot: performance vs sample count and/or event count
- aggregate sortable table
- filters:
  - min events
  - min samples
  - model families
- interpretation text explaining:
  - no universal winner
  - RSF and MultiBranch are strong
  - Cox remains competitive
  - DeepSurv wins only selected cohorts

## 7. Methods & Caveats
Must include:
- concise methodological summary
- cohort omission note for `LAML`
- note that C-index measures discrimination not calibration
- note about small test sets / low event counts
- note on tumor-only modeling vs tumor-normal molecular exploration
- links/text from:
  - `~/tcga_survival/README.md`
  - `~/tcga_survival/METHODS.md`
  - `~/tcga_survival/RESULTS_SUMMARY.md`

# Plot recommendations by data type

Use these plot types:

Cohort metadata:
- bar plots
- stacked bars
- dot plots
- tables

Survival:
- Kaplan–Meier curves
- bar/dot plots for C-index
- heatmaps for cohort × model performance
- scatter for train/test or performance vs event count

RNA / miRNA / methylation / CNV single-feature:
- boxplots
- violin plots
- jitter/beeswarm

Many features:
- heatmaps
- volcano/ranking plots if computed
- PCA / UMAP

Mutation:
- frequency bar plots
- oncoprint-style matrix if feasible
- not boxplots

Avoid:
- pie charts
- overly crowded generic dashboards

# UI/layout expectations

Use `bslib` with a custom theme and add a CSS file in `www/styles.css`.

Layout should be:
- top navbar or sidebar-driven major page navigation
- each major page should feel editorial/scientific
- left controls
- central figure area
- right contextual interpretation/caveat panel when useful

Use:
- `page_navbar()` or similarly polished structure
- nested cards/panels sparingly and elegantly
- avoid clutter

The app should look like a computational biology product, not a business dashboard.

# Technical implementation details

## Data loading
Implement robust path-based loaders in `R/utils_io.R` and `R/utils_data_prep.R`.

Need helpers to:
- discover completed cohort directories
- read `summary.csv` across cohorts
- read aggregate summary
- read cohort metadata from config / run summaries
- locate saved figure files
- read parquet matrices with `arrow`
- parse survival data
- use cache files when available
- attempt tumor/normal extraction from raw source files only when needed

## Caching
Because some raw data are large:
- use `shiny_app/data_cache/` if present
- lazy-load cohort/modality-specific data
- do not load all raw matrices at startup

## Error handling
- missing files should not crash the app
- display informative empty-state panels
- explicitly mark omitted/incomplete cohorts

## Reactivity
- central shared reactive state for selected cohort/modality/grouping where useful
- modular server functions per page
- avoid duplicated data loading

# Success criteria

A successful build should satisfy all of the following:
- app launches locally with `shiny::runApp('~/tcga_survival/shiny_app')`
- Overview page renders with real project data
- Cross-cohort benchmark heatmap renders
- cohort table renders
- at least one cohort works end-to-end on Tumor–Normal Biology page
- at least one clustering workflow works end-to-end on Molecular Structure page
- at least one cohort shows survival-model results and saved figures correctly
- `TCGA-LAML` is visibly marked as excluded with reason
- missing cohorts/files do not crash the app

# Performance requirements

- summary pages should load in roughly < 5–10 sec locally
- cohort-specific pages should lazy-load data
- clustering computations should be limited to selected cohort/modality/feature subset
- expensive computations should be cached or precomputed where practical

# Export/download requirements

Include user-facing export where practical:
- download current summary table as CSV
- download benchmark table as CSV
- download cluster assignments as CSV if clustering is run
- download current plot as PNG if reasonably easy

# Feature annotation requirements

Where possible, support feature metadata columns such as:
- feature id
- modality
- rank / statistic
- tumor-normal summary values

Do not overpromise external annotation if not available, but structure code so it can be added later.

# Reproducibility / session requirements

Where feasible:
- show current analysis settings on clustering pages
- structure code so bookmarking/session state can be added later
- keep cluster method, k, modality, feature subset visible to the user

# Deliverables

1. Full Shiny app scaffold in `~/tcga_survival/shiny_app/`
2. Modularized code, not one giant `app.R`
3. Working README explaining:
   - how to prepare data
   - how to run app
   - package requirements
4. If needed, data prep script(s)
5. Reasonable CSS styling
6. App should run locally with:
   - `Rscript shiny_app/scripts/prepare_app_data.R` if needed
   - then `shiny::runApp('~/tcga_survival/shiny_app')`

# Important content to encode in the app

Scientific conclusions already established:
- completed cohorts: 16
- excluded: `TCGA-LAML`
- no universal best model
- win counts:
  - RSF: 5
  - MultiBranch: 5
  - Cox elastic-net: 4
  - DeepSurv: 2
- multimodal deep learning helps in selected cancers
- classical baselines remain highly competitive

Also surface this cohort-best table from the current aggregate summary:
- ACC: Cox EN
- BLCA: RSF
- BRCA: Cox EN
- CESC: RSF
- CHOL: RSF
- COAD: RSF
- DLBC: MultiBranch
- ESCA: MultiBranch
- GBM: DeepSurv
- HNSC: MultiBranch
- KICH: MultiBranch
- KIRC: Cox EN
- KIRP: MultiBranch
- LGG: DeepSurv
- LIHC: RSF
- LUAD: MultiBranch

# Build priorities

Priority order:
1. app structure + benchmark pages
2. cohort atlas
3. molecular tumor-normal exploration
4. clustering/molecular structure
5. polish and styling

Start by making a clean v1 that works with existing files. Prefer correctness and clarity over too many flashy features.

When done:
- summarize created files
- explain any assumptions
- note any missing data limitations
- suggest next improvements
