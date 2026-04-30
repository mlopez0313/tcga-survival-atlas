Prepare data for an R Shiny app for the project in `~/tcga_survival`.

Important: **do not rebuild the ML preprocessing pipeline from scratch**. Reuse the existing processed outputs and only create additional app-specific caches where necessary.

# Context

Project roots:
- repo: `~/tcga_survival`
- main outputs: `~/mnt/datapool/tcga_survival`
- pan-cohort outputs: `~/mnt/datapool/tcga_survival_cohorts`

The ML pipeline has already produced substantial prepared data for modeling, including:
- per-cohort processed matrices in `TCGA-*/data/processed/*.parquet`
- survival labels in `y_survival.parquet`
- model summaries in `results/metrics/summary.csv`
- saved figures and metrics in `results/figures/` and `results/metrics/`
- aggregate summary files in `_summary/`

This is enough for benchmark/survival pages.

However, the app also needs:
- tumor vs normal molecular exploration
- sample-level clustering
- interrogation of survival and clinical labels by cluster
- flexible molecular views that are likely not in exactly the right app-ready format

# Goal

Create a lightweight app-specific preparation layer that:
- reuses existing processed outputs first
- only derives additional caches where needed
- preserves sample-level tumor/normal structure where possible
- produces fast-to-load files for Shiny

Write outputs under:
- `~/tcga_survival/shiny_app/data_cache/`

If needed, create script:
- `~/tcga_survival/shiny_app/scripts/prepare_app_data.R`

# Principles

- Do not rerun model training
- Do not duplicate huge raw data unnecessarily
- Keep caches compact and app-oriented
- Prefer derived summaries over fully exploded huge tables when possible
- Document all assumptions

# What to prepare

## 1. Cohort metadata cache
Create a compact cohort metadata table, e.g.:
- cohort
- completed status
- excluded status
- exclusion reason
- n_train
- n_test
- events_train
- events_test
- best_model
- best_test_cindex
- modalities available
- paths to figures/results

Use sources like:
- `_summary/aggregate_best_models.csv`
- `TCGA-*/run_summary.json`
- `TCGA-*/results/metrics/summary.csv`
- `TCGA-*/data/processed/`

Save as something like:
- `data_cache/cohort_metadata.rds`
- `data_cache/cohort_metadata.csv`

## 2. Benchmark cache
Create app-friendly cross-cohort benchmark tables:
- all cohort × model metrics
- best-model table
- win counts table
- optional long format for plotting

Save as e.g.:
- `data_cache/benchmark_long.rds`
- `data_cache/benchmark_best.rds`
- `data_cache/win_counts.rds`

## 3. Figure manifest cache
Create a manifest of saved figures for each cohort:
- cohort
- figure path
- figure type
- model

Detect figure types from filenames such as:
- `cindex_comparison.png`
- `km_*.png`
- `deepsurv_training.png`
- `multibranch_training.png`
- `cox_top_coefficients.png`
- `rsf_top_features.png`

Save as:
- `data_cache/figure_manifest.rds`

## 4. Survival modeling cache
Create app-friendly per-cohort metric tables and survival metadata references.
Use existing per-cohort `summary.csv` and `y_survival.parquet`.

Save compact objects needed for fast loading in app.

## 5. Sample annotation cache for molecular exploration
This is critical.
Create a sample-level annotation table per cohort where possible, including:
- sample_id
- patient_id if derivable
- cohort
- sample_type / tumor_normal status
- survival time/event if matched
- stage
- sex
- age if available
- modality availability flags

You may need to derive this from raw files and/or processed data plus clinical tables.

This does not need to be perfect for all cohorts, but should be robust and well documented.

Save as:
- `data_cache/sample_annotations.rds`
- or one file per cohort if easier

## 6. Modality-specific app caches
For interactive tumor-vs-normal plotting and clustering, prepare manageable sample-level matrices or summaries for:
- RNA
- miRNA
- methylation
- CNV
- mutation if feasible

Important:
- do not cache the entire raw universe in absurdly large format if unnecessary
- prefer top variable / most useful features for app interactivity
- keep enough information for:
  - box/violin plots
  - PCA/UMAP
  - clustering
  - cluster-feature heatmaps

Suggested strategy:
For each cohort × modality, prepare:
- matrix of top N variable features (e.g. 500-2000 depending on modality)
- sample annotation joinable to matrix columns/rows
- feature metadata table
- optional long-format subset for quick feature plots

Examples:
- `data_cache/modalities/TCGA-LUAD_rna_topvar.rds`
- `data_cache/modalities/TCGA-LUAD_methylation_topvar.rds`
- etc.

Potential defaults:
- RNA: top 1000 variable genes
- miRNA: top 200-500
- methylation: top 1000-2000 probes for app exploration
- CNV: top 500-1000 features
- mutation: top mutated genes / binary matrix

If tumor-normal labels are not preserved in processed matrices, use raw files to recover sample-level structure.

## 7. Tumor vs normal feature summary cache
For each cohort/modality, if possible precompute feature-level summary stats comparing tumor vs normal:
- mean_tumor
- mean_normal
- median_tumor
- median_normal
- difference
- logFC-like measure where appropriate
- p-value optional if easy
- missingness
- n_tumor
- n_normal

This will make app feature ranking much faster.

Save as e.g.:
- `data_cache/feature_summaries/<cohort>_<modality>_tumor_normal.rds`

## 8. Clustering-ready cache
Prepare app-ready clustering inputs per cohort/modality:
- top-variable feature matrices
- scaled matrices if useful
- optional precomputed PCA/UMAP

At minimum the app should support clustering methods like:
- hierarchical clustering
- k-means

Optional but nice:
- precompute PCA coordinates
- precompute UMAP for default settings

Save e.g.:
- `data_cache/clustering/<cohort>_<modality>_matrix.rds`
- `data_cache/clustering/<cohort>_<modality>_pca.rds`
- `data_cache/clustering/<cohort>_<modality>_umap.rds`

# Special notes

## LAML
- mark `TCGA-LAML` as excluded
- include exclusion reason in metadata cache
- do not fail the prep because of LAML

## Existing processed data reuse
Use the existing processed outputs for:
- survival benchmark views
- per-cohort modeling summaries
- train/test sizes and events

Only derive extra data where the app needs sample-level or tumor/normal exploration that the modeling pipeline did not preserve conveniently.

# Package suggestions
Use R packages such as:
- `data.table`
- `dplyr`
- `tidyr`
- `purrr`
- `stringr`
- `jsonlite`
- `readr`
- `arrow`
- `fs`
- `matrixStats`
- `uwot` if precomputing UMAP
- `stats`

# Success criteria

A successful prep pass should:
- create `shiny_app/data_cache/`
- write cache files for benchmark pages without rerunning ML
- create at least one working cohort/modality cache for tumor-normal exploration
- create at least one working clustering-ready cache
- include `TCGA-LAML` as excluded with reason
- not crash if some cohorts/modalities are incomplete
- print a summary of generated artifacts

# Performance requirements

- do not materialize giant unnecessary long tables for all raw features if avoidable
- prefer lazy cohort/modality-specific cache files
- keep first app load feasible (< ~5–10 sec locally for summary pages)
- cache expensive derived computations

# Deliverables

1. Create `shiny_app/scripts/prepare_app_data.R`
2. Write cache outputs under `shiny_app/data_cache/`
3. Create a short manifest/README describing each cache file and how it was derived
4. Print a summary of:
   - cohorts discovered
   - completed vs excluded
   - caches written
   - assumptions / limitations

# Important implementation detail

Please be conservative and pragmatic:
- prefer a working, fast v1 cache layer
- do not overengineer
- do not create huge unnecessarily denormalized files
- if some cohorts/modalities are too awkward for tumor-normal exploration in v1, degrade gracefully and document it

At the end, summarize:
- what cache files were created
- which app pages they support
- what limitations remain
