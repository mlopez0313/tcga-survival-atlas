# utils_config.R --------------------------------------------------------------
# Single source of truth for filesystem paths and runtime configuration.
#
# All paths are resolved from env vars first, then sensible defaults, so the
# same code can run on a developer laptop, a Posit Connect host, or inside a
# container without touching the modules.
#
# Recognised environment variables:
#   TCGA_SURVIVAL_REPO_ROOT      Repo root (defaults to parent of shiny_app/).
#   TCGA_SURVIVAL_APP_ROOT       App directory (defaults to shiny_app/).
#   TCGA_SURVIVAL_CACHE_DIR      data_cache/ override.
#   TCGA_SURVIVAL_COHORT_ROOT    Per-cohort raw + processed outputs root.
#   TCGA_SURVIVAL_SUMMARY_DIR    Aggregate summary directory.
#   TCGA_SURVIVAL_LOG_LEVEL      info | warn | error (default: info).
# -----------------------------------------------------------------------------

#' Tiny logger that prefixes every line so it shows up clearly in the
#' Shiny / R console without any external dependency.
log_msg <- function(..., level = c("info","warn","error","debug")) {
  level <- match.arg(level)
  threshold <- tolower(Sys.getenv("TCGA_SURVIVAL_LOG_LEVEL", "info"))
  ranks <- c(debug = 1L, info = 2L, warn = 3L, error = 4L)
  if (ranks[[level]] < ranks[[threshold %||% "info"]]) return(invisible())
  tag <- sprintf("[%s][tcga-atlas] ",
                 format(Sys.time(), "%H:%M:%S"))
  msg <- paste0(tag, paste0(..., collapse = ""))
  if (level == "error") message(msg)
  else if (level == "warn") warning(msg, call. = FALSE)
  else message(msg)
  invisible()
}

# `%||%` is defined in utils_data_prep.R; redefine defensively here so this
# file is self-sufficient even if sourced first.
if (!exists("%||%", inherits = TRUE)) {
  `%||%` <- function(a, b) {
    if (length(a) == 0L) return(b)
    if (length(a) == 1L) {
      if (is.na(a)) return(b)
      if (is.character(a) && !nzchar(a)) return(b)
    }
    a
  }
}

#' Resolve a path from an env var or fall back to a default. Optional
#' `must_exist` only emits a warning, never a hard error — the goal is for
#' the app to keep running even when peripheral dirs are missing.
.resolve_path <- function(env_var, default, must_exist = FALSE,
                           label = env_var) {
  raw <- Sys.getenv(env_var, unset = NA_character_)
  p   <- if (!is.na(raw) && nzchar(raw)) raw else default
  p   <- normalizePath(p, winslash = "/", mustWork = FALSE)
  if (must_exist && !file.exists(p)) {
    log_msg(sprintf("path missing: %s = %s (override with $%s)",
                    label, p, env_var), level = "warn")
  }
  p
}

#' Build the canonical APP_PATHS object. Called once from global.R; the
#' returned list is also re-exported as legacy all-caps constants for
#' backwards compatibility with existing module code.
build_app_paths <- function(app_root_hint = NULL) {
  # 1. anchor: app_root --------------------------------------------------
  default_app_root <- if (!is.null(app_root_hint) &&
                          basename(app_root_hint) == "shiny_app")
    app_root_hint
  else if (basename(getwd()) == "shiny_app") getwd()
  else file.path(Sys.getenv("HOME"), "tcga_survival", "shiny_app")
  app_root  <- .resolve_path("TCGA_SURVIVAL_APP_ROOT", default_app_root,
                              must_exist = TRUE, label = "app_root")

  # 2. derived defaults --------------------------------------------------
  default_repo_root    <- normalizePath(file.path(app_root, ".."),
                                          winslash = "/", mustWork = FALSE)
  default_cache_dir    <- file.path(app_root, "data_cache")
  default_cohort_root  <- file.path(Sys.getenv("HOME"),
                                     "mnt/datapool/tcga_survival_cohorts")
  default_summary_dir  <- file.path(default_cohort_root, "_summary")

  repo_root   <- .resolve_path("TCGA_SURVIVAL_REPO_ROOT",   default_repo_root)
  cache_dir   <- .resolve_path("TCGA_SURVIVAL_CACHE_DIR",   default_cache_dir,
                                 must_exist = TRUE, label = "cache_dir")
  cohort_root <- .resolve_path("TCGA_SURVIVAL_COHORT_ROOT", default_cohort_root,
                                 must_exist = FALSE, label = "cohort_root")
  summary_dir <- .resolve_path("TCGA_SURVIVAL_SUMMARY_DIR", default_summary_dir,
                                 must_exist = FALSE, label = "summary_dir")

  list(
    app_root       = app_root,
    repo_root      = repo_root,
    cache_dir      = cache_dir,
    cohort_root    = cohort_root,
    summary_dir    = summary_dir,
    sample_anns    = file.path(cache_dir, "sample_annotations"),
    modalities     = file.path(cache_dir, "modalities"),
    feat_summ      = file.path(cache_dir, "feature_summaries"),
    clustering     = file.path(cache_dir, "clustering"),
    cohort_summ    = file.path(cache_dir, "cohort_summaries"),
    manifest_md    = file.path(cache_dir, "MANIFEST.md")
  )
}

#' Run the full battery of startup checks. Returns a data.table of
#' (check, status, detail). Status is one of "ok"/"warn"/"missing".
#' Never throws — it is safe to call from anywhere.
validate_app_data <- function(paths = APP_PATHS, app_data = NULL) {
  rows <- list()
  add <- function(check, status, detail = "") {
    rows[[length(rows) + 1L]] <<- data.table::data.table(
      check = check, status = status, detail = detail)
  }

  # 1. directories ------------------------------------------------------
  for (nm in c("app_root","cache_dir","cohort_root","summary_dir",
               "sample_anns","modalities","feat_summ","clustering",
               "cohort_summ")) {
    p <- paths[[nm]]
    if (is.null(p) || !nzchar(p)) { add(nm, "missing", "(unset)"); next }
    if (!file.exists(p))           { add(nm, "missing", p);        next }
    add(nm, "ok", p)
  }

  # 2. headline cache files --------------------------------------------
  needed <- list(
    cohort_metadata  = "cohort_metadata.rds",
    benchmark_long   = "benchmark_long.rds",
    benchmark_best   = "benchmark_best.rds",
    win_counts       = "win_counts.rds",
    figure_manifest  = "figure_manifest.rds",
    survival_summary = "survival_summary.rds",
    sample_annotations = "sample_annotations.rds",
    manifest         = "MANIFEST.md"
  )
  for (k in names(needed)) {
    p <- file.path(paths$cache_dir, needed[[k]])
    if (file.exists(p))
      add(paste0("cache:", k), "ok", basename(p))
    else
      add(paste0("cache:", k), "missing",
          sprintf("missing %s — run scripts/prepare_app_data.R",
                  basename(p)))
  }

  # 3. cohort discovery ------------------------------------------------
  meta <- if (!is.null(app_data) && !is.null(app_data$cohort_metadata))
            app_data$cohort_metadata
          else tryCatch(readRDS(file.path(paths$cache_dir,
                                           "cohort_metadata.rds")),
                        error = function(e) NULL)
  if (is.null(meta)) {
    add("cohorts:metadata", "missing", "no cohort_metadata.rds")
  } else {
    n_total <- nrow(meta)
    n_done  <- sum(meta$status == "completed", na.rm = TRUE)
    add("cohorts:metadata", "ok",
        sprintf("%d cohorts (%d completed)", n_total, n_done))
  }

  # 4. modality cache coverage -----------------------------------------
  mod_files <- list.files(paths$modalities, pattern = "\\.rds$")
  add("cohorts:modality_cache",
      if (length(mod_files) > 0) "ok" else "missing",
      sprintf("%d modality cache files", length(mod_files)))

  cl_files <- list.files(paths$clustering, pattern = "\\.rds$")
  add("cohorts:clustering_cache",
      if (length(cl_files) > 0) "ok" else "warn",
      sprintf("%d PCA/UMAP cache files", length(cl_files)))

  # 5. figure paths (sample 5 cohorts only — cheap O(1) per check) -----
  if (!is.null(meta)) {
    sample_cohorts <- head(as.character(meta[status == "completed", cohort]), 5L)
    for (co in sample_cohorts) {
      figs <- list.files(file.path(paths$cohort_root, co, "results", "figures"),
                          pattern = "\\.png$", full.names = TRUE)
      add(paste0("figures:", co),
          if (length(figs) > 0) "ok" else "warn",
          sprintf("%d png file%s", length(figs),
                   if (length(figs) == 1L) "" else "s"))
    }
  }

  # 6. one-cohort end-to-end smoke test -------------------------------
  smoke_co <- if (!is.null(meta) && "TCGA-LUAD" %in% meta$cohort) "TCGA-LUAD"
              else if (!is.null(meta)) {
                m <- meta[status == "completed"]
                if (nrow(m)) as.character(m$cohort[1]) else NA_character_
              } else NA_character_
  if (!is.na(smoke_co)) {
    sa <- tryCatch(readRDS(file.path(paths$sample_anns,
                                       paste0(smoke_co, ".rds"))),
                    error = function(e) NULL)
    add(paste0("smoke:", smoke_co, ":sample_annotation"),
        if (!is.null(sa) && nrow(sa)) "ok" else "missing",
        if (!is.null(sa)) sprintf("%d samples", nrow(sa)) else "no annotation")

    mod_pat <- sprintf("^%s_(.+)_topvar\\.rds$", smoke_co)
    mods    <- sub(mod_pat, "\\1", list.files(paths$modalities, pattern = mod_pat))
    add(paste0("smoke:", smoke_co, ":modalities"),
        if (length(mods) > 0) "ok" else "missing",
        if (length(mods)) paste(mods, collapse = ", ") else "no modality cache")

    if (length(mods)) {
      mo_path <- file.path(paths$modalities,
                            sprintf("%s_%s_topvar.rds", smoke_co, mods[1]))
      mo <- tryCatch(readRDS(mo_path), error = function(e) NULL)
      ok <- !is.null(mo) && !is.null(mo$matrix) && ncol(mo$matrix) >= 4L
      add(paste0("smoke:", smoke_co, ":clusterable"),
          if (ok) "ok" else "warn",
          if (ok) sprintf("%s: %d features × %d samples",
                          mods[1], nrow(mo$matrix), ncol(mo$matrix))
          else "matrix unavailable or too small to cluster")
    }
  }

  out <- data.table::rbindlist(rows, fill = TRUE)
  attr(out, "n_warn")    <- sum(out$status == "warn",    na.rm = TRUE)
  attr(out, "n_missing") <- sum(out$status == "missing", na.rm = TRUE)
  out
}

#' Pretty-print validation results to the console with a one-line summary.
log_validation <- function(v) {
  if (is.null(v) || !nrow(v)) {
    log_msg("validate_app_data() returned no rows", level = "warn"); return(invisible())
  }
  n_ok      <- sum(v$status == "ok",      na.rm = TRUE)
  n_warn    <- sum(v$status == "warn",    na.rm = TRUE)
  n_missing <- sum(v$status == "missing", na.rm = TRUE)
  log_msg(sprintf("validation: %d ok · %d warn · %d missing",
                  n_ok, n_warn, n_missing))
  for (i in which(v$status %in% c("warn","missing"))) {
    log_msg(sprintf("  %-7s %-32s %s",
                    v$status[i], v$check[i], v$detail[i]),
            level = "warn")
  }
  invisible()
}
