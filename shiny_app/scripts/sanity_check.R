#!/usr/bin/env Rscript
# scripts/sanity_check.R ------------------------------------------------------
# Headless sanity check for the TCGA Survival Atlas app.
#
# Runs the same `validate_app_data()` battery the Methods → Diagnostics card
# uses, plus a representative cohort end-to-end exercise (load metadata →
# load one modality cache → cluster the matrix → compute interpretation).
#
# Usage:
#   Rscript shiny_app/scripts/sanity_check.R
#   TCGA_SURVIVAL_CACHE_DIR=/srv/data Rscript shiny_app/scripts/sanity_check.R
#
# Exit status:
#   0 — all checks ok or merely "warn"
#   1 — at least one "missing" check
#   2 — unhandled R error
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L) return(b)
  if (length(a) == 1L && (is.na(a) || (is.character(a) && !nzchar(a))))
    return(b)
  a
}

# Locate the app dir relative to this script ---------------------------------
script_arg <- (function() {
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(ofile)) return(ofile)
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) return(sub("^--file=", "", m[1]))
  "shiny_app/scripts/sanity_check.R"
})()
this_file <- normalizePath(script_arg, mustWork = FALSE)
app_dir <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
if (!file.exists(file.path(app_dir, "global.R"))) {
  app_dir <- normalizePath(file.path(Sys.getenv("HOME"),
                                       "tcga_survival", "shiny_app"),
                            mustWork = FALSE)
}
if (!file.exists(file.path(app_dir, "global.R"))) {
  stop("Could not locate shiny_app/global.R — pass the right working dir.")
}

setwd(app_dir)
cat("[sanity] app dir: ", app_dir, "\n", sep = "")

# Source utilities only (don't start the Shiny app) --------------------------
source("R/utils_config.R", local = FALSE)
source("R/utils_data_prep.R", local = FALSE)
source("R/utils_io.R", local = FALSE)

paths    <- build_app_paths(app_dir)
app_data <- load_app_data(paths$cache_dir)

cat(sprintf("[sanity] app_root    = %s\n", paths$app_root))
cat(sprintf("[sanity] cache_dir   = %s\n", paths$cache_dir))
cat(sprintf("[sanity] cohort_root = %s\n", paths$cohort_root))
cat(sprintf("[sanity] summary_dir = %s\n", paths$summary_dir))

v <- tryCatch(validate_app_data(paths, app_data),
              error = function(e) {
                cat("[sanity] FATAL: ", conditionMessage(e), "\n")
                quit(status = 2)
              })

# Pretty-print the validation result ----------------------------------------
fmt_row <- function(check, status, detail) {
  pad <- function(s, n) formatC(s, width = -n, flag = "-")
  status <- pad(status, 7)
  check  <- pad(check, 36)
  sprintf("  %s  %s  %s", status, check, detail)
}
cat("[sanity] validation rows:\n")
for (i in seq_len(nrow(v))) {
  cat(fmt_row(v$check[i], v$status[i], v$detail[i]), "\n")
}

n_ok      <- sum(v$status == "ok",      na.rm = TRUE)
n_warn    <- sum(v$status == "warn",    na.rm = TRUE)
n_missing <- sum(v$status == "missing", na.rm = TRUE)
cat(sprintf("\n[sanity] summary: %d ok · %d warn · %d missing\n\n",
            n_ok, n_warn, n_missing))

# Representative-cohort end-to-end probe -------------------------------------
meta <- app_data$cohort_metadata
if (!is.null(meta)) {
  smoke <- if ("TCGA-LUAD" %in% meta$cohort) "TCGA-LUAD"
           else as.character(meta[status == "completed", cohort][1])
  if (!is.na(smoke)) {
    cat(sprintf("[sanity] e2e probe: %s\n", smoke))
    sa <- load_sample_annotation(smoke, app_data)
    cat(sprintf("  sample_annotation rows  : %s\n",
                if (is.null(sa)) "MISSING" else nrow(sa)))
    mods <- modalities_with_cache(smoke, app_data)
    cat(sprintf("  modality caches present : %s\n",
                if (length(mods)) paste(mods, collapse = ", ") else "NONE"))
    if (length(mods)) {
      mo <- load_modality_cache(smoke, mods[1], app_data)
      ok <- !is.null(mo) && !is.null(mo$matrix)
      cat(sprintf("  primary matrix          : %s\n",
                  if (ok) sprintf("%d × %d (%s)",
                                  nrow(mo$matrix), ncol(mo$matrix), mods[1])
                  else "MISSING"))
      if (ok && ncol(mo$matrix) >= 6L) {
        cl <- tryCatch(run_clustering(mo$matrix, method = "kmeans", k = 3L,
                                       scale_features = TRUE),
                       error = function(e) NULL)
        cat(sprintf("  k-means k=3 clusters    : %s\n",
                    if (!is.null(cl)) paste(table(cl), collapse = "/")
                    else "FAILED"))
      }
    }
    tn <- load_tn_summary(smoke, mods[1] %||% "rna", app_data)
    cat(sprintf("  tumor/normal summary    : %s\n",
                if (is.null(tn)) "n/a (no normals or modality)"
                else sprintf("%d features", nrow(tn))))
  }
}

quit(status = if (n_missing > 0L) 1L else 0L)
