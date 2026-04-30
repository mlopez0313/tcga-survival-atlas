# utils_io.R ------------------------------------------------------------------
# Cache + raw IO helpers. All functions here are pure (no Shiny reactivity).
# -----------------------------------------------------------------------------

#' Safe RDS reader that returns NULL on missing/error.
read_rds_safe <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

#' Load all light caches. Returns a list of data.tables (or NULL where missing).
load_app_data <- function(cache_dir) {
  list(
    cache_dir         = cache_dir,
    cohort_metadata   = read_rds_safe(file.path(cache_dir, "cohort_metadata.rds")),
    benchmark_long    = read_rds_safe(file.path(cache_dir, "benchmark_long.rds")),
    benchmark_best    = read_rds_safe(file.path(cache_dir, "benchmark_best.rds")),
    win_counts        = read_rds_safe(file.path(cache_dir, "win_counts.rds")),
    figure_manifest   = read_rds_safe(file.path(cache_dir, "figure_manifest.rds")),
    survival_summary  = read_rds_safe(file.path(cache_dir, "survival_summary.rds")),
    sample_anns_path  = file.path(cache_dir, "sample_annotations"),
    modalities_path   = file.path(cache_dir, "modalities"),
    feat_summ_path    = file.path(cache_dir, "feature_summaries"),
    clustering_path   = file.path(cache_dir, "clustering"),
    cohort_summ_path  = file.path(cache_dir, "cohort_summaries"),
    manifest_md       = file.path(cache_dir, "MANIFEST.md")
  )
}

#' Lazy: per-cohort sample annotation.
load_sample_annotation <- function(cohort, app_data = APP_DATA) {
  read_rds_safe(file.path(app_data$sample_anns_path, paste0(cohort, ".rds")))
}

#' Lazy: per-cohort y_survival table.
load_cohort_survival <- function(cohort, app_data = APP_DATA) {
  read_rds_safe(file.path(app_data$cohort_summ_path,
                          paste0(cohort, "_survival.rds")))
}

#' Lazy: per-cohort/modality top-variable matrix cache.
#' Returns list(matrix, feature_meta, sample_meta, modality, cohort) or NULL.
load_modality_cache <- function(cohort, modality, app_data = APP_DATA) {
  read_rds_safe(file.path(app_data$modalities_path,
                          sprintf("%s_%s_topvar.rds", cohort, modality)))
}

#' Lazy: per-cohort/modality tumor-vs-normal feature stats.
load_tn_summary <- function(cohort, modality, app_data = APP_DATA) {
  read_rds_safe(file.path(app_data$feat_summ_path,
                          sprintf("%s_%s_tumor_normal.rds", cohort, modality)))
}

#' Lazy: per-cohort/modality precomputed PCA scores.
load_pca_cache <- function(cohort, modality, app_data = APP_DATA) {
  read_rds_safe(file.path(app_data$clustering_path,
                          sprintf("%s_%s_pca.rds", cohort, modality)))
}

#' Lazy: per-cohort/modality precomputed UMAP scores.
load_umap_cache <- function(cohort, modality, app_data = APP_DATA) {
  read_rds_safe(file.path(app_data$clustering_path,
                          sprintf("%s_%s_umap.rds", cohort, modality)))
}

#' Per-cohort metric summary.csv (raw read, lazy). Resolved via APP_PATHS so
#' the cohort root can be overridden through the TCGA_SURVIVAL_COHORT_ROOT env
#' var without touching this file.
load_cohort_summary_csv <- function(cohort, paths = APP_PATHS) {
  p <- file.path(paths$cohort_root, cohort, "results", "metrics", "summary.csv")
  if (!file.exists(p)) return(NULL)
  tryCatch(data.table::fread(p), error = function(e) NULL)
}

#' List saved figures for a cohort (full paths). Returns character(0) when
#' the figures directory does not exist (e.g. excluded cohort, or a freshly
#' bootstrapped deployment without the modeling outputs).
list_cohort_figures <- function(cohort, paths = APP_PATHS) {
  fdir <- file.path(paths$cohort_root, cohort, "results", "figures")
  if (!dir.exists(fdir)) return(character())
  list.files(fdir, pattern = "\\.png$", full.names = TRUE)
}

#' Raw text file reader that returns "" on failure (used for METHODS/README).
read_text_safe <- function(path) {
  if (!file.exists(path)) return("")
  tryCatch(paste(readLines(path, warn = FALSE), collapse = "\n"),
           error = function(e) "")
}

#' Available cohorts based on cohort_metadata; ordered, with status flags.
get_completed_cohorts <- function(app_data = APP_DATA) {
  m <- app_data$cohort_metadata
  if (is.null(m)) return(character())
  as.character(m[status == "completed", cohort])
}

get_all_cohorts <- function(app_data = APP_DATA) {
  m <- app_data$cohort_metadata
  if (is.null(m)) return(character())
  as.character(m[, cohort])
}

#' Modalities present (both processed and via raw data) for a given cohort,
#' restricted to those that actually have a modality cache file written.
modalities_with_cache <- function(cohort, app_data = APP_DATA) {
  pat <- sprintf("^%s_(.+)_topvar\\.rds$", cohort)
  files <- list.files(app_data$modalities_path, pattern = pat)
  if (!length(files)) return(character())
  sub(pat, "\\1", files)
}
