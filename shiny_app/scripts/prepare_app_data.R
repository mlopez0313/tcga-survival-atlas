#!/usr/bin/env Rscript
# =============================================================================
# prepare_app_data.R
#
# Build the lightweight, app-specific cache layer that the TCGA-survival
# Shiny app loads at startup. Reuses the existing modeling outputs first;
# only derives extra caches where the app needs sample-level / tumor-normal
# views the modeling pipeline did not preserve.
#
# Inputs (read-only):
#   ~/mnt/datapool/tcga_survival_cohorts/<TCGA-XXXX>/...
#   ~/mnt/datapool/tcga_survival_cohorts/_summary/aggregate_best_models.csv
#
# Outputs:
#   ~/tcga_survival/shiny_app/data_cache/
#       cohort_metadata.{rds,csv}
#       benchmark_long.rds
#       benchmark_best.rds
#       win_counts.rds
#       figure_manifest.rds
#       survival_summary.rds
#       sample_annotations.rds
#       sample_annotations/<cohort>.rds
#       modalities/<cohort>_<modality>_topvar.rds
#       feature_summaries/<cohort>_<modality>_tumor_normal.rds
#       clustering/<cohort>_<modality>_pca.rds
#       MANIFEST.md
#
# Usage:
#   Rscript shiny_app/scripts/prepare_app_data.R [--quick] [--cohorts=LUAD,LUSC]
#
# Flags:
#   --quick               only build light caches (skip raw modality work)
#   --cohorts=A,B,C       restrict heavy modality work to these cohorts
#                         (cohort metadata / benchmarks always cover all)
#   --no-methylation      skip the (large) methylation raw matrix
#   --no-umap             skip UMAP precompute (PCA still computed)
#   --features-rna=N      override top-variable RNA gene count (default 1000)
#   --features-meth=N     override top-variable methylation probes (default 2000)
#   --features-cnv=N      override top-variable CNV features  (default 1000)
#   --features-mirna=N    override top-variable miRNAs       (default 200)
#   --features-mut=N      override top mutated genes         (default 300)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(stringr)
})

# ------------------------------------------------------------ paths & options
# All paths are overridable via env vars (see shiny_app/R/utils_config.R for
# the shared resolution rules — kept duplicated here so this script remains
# usable without sourcing the app's globals).
.env_or <- function(var, default) {
  v <- Sys.getenv(var, unset = NA_character_)
  if (!is.na(v) && nzchar(v)) v else default
}

HOME           <- Sys.getenv("HOME")
REPO           <- .env_or("TCGA_SURVIVAL_REPO_ROOT",   file.path(HOME, "tcga_survival"))
APP_DIR        <- .env_or("TCGA_SURVIVAL_APP_ROOT",    file.path(REPO, "shiny_app"))
CACHE_DIR      <- .env_or("TCGA_SURVIVAL_CACHE_DIR",   file.path(APP_DIR, "data_cache"))
COHORT_ROOT    <- .env_or("TCGA_SURVIVAL_COHORT_ROOT", file.path(HOME, "mnt/datapool/tcga_survival_cohorts"))
SUMMARY_DIR    <- .env_or("TCGA_SURVIVAL_SUMMARY_DIR", file.path(COHORT_ROOT, "_summary"))
SCRIPT_DIR     <- file.path(APP_DIR, "scripts")
PY_BIN         <- .env_or("TCGA_SURVIVAL_PYTHON",      file.path(HOME, "venv/bin/python"))
PY_HELPER      <- file.path(SCRIPT_DIR, "_pq2csv.py")

EXCLUDED_COHORTS <- list(
  "TCGA-LAML" = paste0(
    "TCGAbiolinks downloader could not retrieve required clinical fields ",
    "(notably `disease_response`); cohort never reached preprocessing."
  )
)

# ----------------------------------------------------- CLI flags --------------
args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(name, default = FALSE) {
  any(args == paste0("--", name)) || default
}
get_opt <- function(name, default = NULL, as_int = FALSE) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  v <- sub(paste0("^--", name, "="), "", hit[1])
  if (as_int) as.integer(v) else v
}

OPT_QUICK         <- get_flag("quick")
OPT_NO_METH       <- get_flag("no-methylation")
OPT_NO_UMAP       <- get_flag("no-umap")
OPT_COHORTS       <- get_opt("cohorts", default = NULL)
N_RNA             <- get_opt("features-rna",   1000L, as_int = TRUE)
N_METH            <- get_opt("features-meth",  2000L, as_int = TRUE)
N_CNV             <- get_opt("features-cnv",   1000L, as_int = TRUE)
N_MIRNA           <- get_opt("features-mirna",  200L, as_int = TRUE)
N_MUT             <- get_opt("features-mut",    300L, as_int = TRUE)

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
for (sub in c("modalities", "feature_summaries", "clustering",
              "sample_annotations", "cohort_summaries")) {
  dir.create(file.path(CACHE_DIR, sub), showWarnings = FALSE, recursive = TRUE)
}

logmsg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), ..., "\n", sep = " ")

# ----------------------------------------------------- parquet helper ---------
# Convert a parquet file to a data.table by piping through the python helper.
read_parquet_dt <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (!file.exists(PY_BIN) || !file.exists(PY_HELPER)) {
    stop("python parquet helper not available: ", PY_BIN, " / ", PY_HELPER)
  }
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)
  status <- system2(PY_BIN, args = c(shQuote(PY_HELPER), shQuote(path), shQuote(tmp)),
                    stdout = FALSE, stderr = FALSE)
  if (status != 0 || !file.exists(tmp)) {
    warning("parquet read failed for ", path)
    return(NULL)
  }
  data.table::fread(tmp)
}

# ----------------------------------------------------- TCGA helpers -----------
tcga_patient_from_sample <- function(s) {
  # TCGA-XX-YYYY-01A-... -> TCGA-XX-YYYY
  v <- strsplit(s, "-", fixed = TRUE)
  vapply(v, function(p) if (length(p) >= 3) paste(p[1:3], collapse = "-") else NA_character_,
         character(1))
}

tcga_sample_type <- function(s) {
  # Sample-type code = 4th block (e.g. "01A"). 01-09 = tumor, 10-19 = normal.
  v <- strsplit(s, "-", fixed = TRUE)
  codes <- vapply(v, function(p) if (length(p) >= 4) substr(p[4], 1, 2) else NA_character_,
                  character(1))
  out <- rep(NA_character_, length(codes))
  num <- suppressWarnings(as.integer(codes))
  out[!is.na(num) & num >= 1 & num <= 9]   <- "tumor"
  out[!is.na(num) & num >= 10 & num <= 19] <- "normal"
  out[!is.na(num) & num >= 20]             <- "control"
  out
}

# Truncate full aliquot barcode to vial-level "sample_short" id:
#   TCGA-XX-YYYY-01A-...  ->  TCGA-XX-YYYY-01A
# Different aliquots from the same vial collapse to one short id; this is
# the level at which we treat samples as "the same biological specimen"
# across modalities.
tcga_sample_short <- function(s) {
  v <- strsplit(s, "-", fixed = TRUE)
  vapply(v, function(p) if (length(p) >= 4) paste(p[1:4], collapse = "-") else NA_character_,
         character(1))
}

# ----------------------------------------------------- discovery --------------
discover_cohorts <- function() {
  dirs <- list.dirs(COHORT_ROOT, recursive = FALSE)
  dirs <- dirs[grepl("/TCGA-[A-Z]+$", dirs)]
  data.table(
    cohort   = basename(dirs),
    root_dir = dirs
  )[order(cohort)]
}

cohort_paths <- function(cohort) {
  base <- file.path(COHORT_ROOT, cohort)
  list(
    cohort       = cohort,
    root         = base,
    raw_dir      = file.path(base, "data/raw"),
    processed    = file.path(base, "data/processed"),
    summary_csv  = file.path(base, "results/metrics/summary.csv"),
    figures_dir  = file.path(base, "results/figures"),
    metrics_dir  = file.path(base, "results/metrics"),
    run_summary  = file.path(base, "run_summary.json"),
    config_yaml  = file.path(base, "config.yaml")
  )
}

cohort_processed_files <- function(cohort) {
  p <- cohort_paths(cohort)$processed
  files <- list.files(p, pattern = "^X_.*\\.parquet$", full.names = TRUE)
  basename(files)
}

cohort_modalities_present <- function(cohort) {
  files <- cohort_processed_files(cohort)
  modalities <- c(
    rna         = "X_rna.parquet",
    mirna       = "X_mirna.parquet",
    methylation = "X_methylation.parquet",
    cnv         = "X_cnv.parquet",
    mutation    = "X_mutation.parquet",
    clinical    = "X_clinical.parquet"
  )
  setNames(modalities %in% files, names(modalities))
}

cohort_raw_files <- function(cohort) {
  raw <- cohort_paths(cohort)$raw_dir
  if (!dir.exists(raw)) return(list())
  files <- list.files(raw, pattern = paste0("^", cohort, "\\."), full.names = TRUE)
  ext_map <- list(
    clinical    = paste0(cohort, ".clinical.tsv.gz"),
    rnaseq_cnt  = paste0(cohort, ".htseq_counts.tsv.gz"),
    rnaseq_fpkm = paste0(cohort, ".htseq_fpkm.tsv.gz"),
    methylation = paste0(cohort, ".methylation450.tsv.gz"),
    mirna       = paste0(cohort, ".mirna.tsv.gz"),
    cnv         = paste0(cohort, ".gistic.tsv.gz"),
    mutation    = paste0(cohort, ".mutect2_snv.tsv.gz")
  )
  out <- lapply(ext_map, function(name) {
    p <- file.path(raw, name)
    if (file.exists(p)) p else NA_character_
  })
  out
}

# =============================================================================
# 1) Cohort metadata + benchmark + figure manifest + survival caches
# =============================================================================
build_cohort_metadata <- function(cohorts) {
  agg_csv <- file.path(SUMMARY_DIR, "aggregate_best_models.csv")
  agg <- if (file.exists(agg_csv)) data.table::fread(agg_csv) else
    data.table(cohort = character(), best_model = character(),
               test_c_index = numeric(), n_train = integer(),
               n_test = integer(), events_train = integer(),
               events_test = integer())

  rows <- lapply(cohorts$cohort, function(co) {
    p <- cohort_paths(co)
    excluded   <- co %in% names(EXCLUDED_COHORTS)
    excl_rsn   <- if (excluded) EXCLUDED_COHORTS[[co]] else NA_character_
    has_proc   <- file.exists(file.path(p$processed, "y_survival.parquet"))
    has_summ   <- file.exists(p$summary_csv)
    has_figs   <- dir.exists(p$figures_dir) &&
                  length(list.files(p$figures_dir, pattern = "\\.png$")) > 0

    completed  <- has_proc && has_summ
    status     <- if (excluded)        "excluded"
                  else if (completed)  "completed"
                  else if (has_proc)   "partial"
                  else                 "incomplete"

    run <- list()
    if (file.exists(p$run_summary)) {
      run <- tryCatch(jsonlite::fromJSON(p$run_summary), error = function(e) list())
    }

    # n_train/test/events: prefer summary.csv (has it per model, all the same)
    n_train <- NA_integer_; n_test <- NA_integer_
    e_train <- NA_integer_; e_test <- NA_integer_
    best_model <- NA_character_; best_c <- NA_real_
    if (has_summ) {
      s <- tryCatch(data.table::fread(p$summary_csv), error = function(e) NULL)
      if (!is.null(s) && nrow(s) > 0) {
        n_train <- as.integer(s$n_train[1])
        n_test  <- as.integer(s$n_test[1])
        e_train <- as.integer(s$events_train[1])
        e_test  <- as.integer(s$events_test[1])
      }
    }
    a <- agg[cohort == co]
    if (nrow(a) > 0) {
      best_model <- a$best_model[1]
      best_c     <- a$test_c_index[1]
    }

    mods <- cohort_modalities_present(co)
    list(
      cohort                = co,
      status                = status,
      excluded              = excluded,
      exclusion_reason      = excl_rsn,
      n_train               = n_train,
      n_test                = n_test,
      events_train          = e_train,
      events_test           = e_test,
      best_model            = best_model,
      best_test_cindex      = best_c,
      modality_rna          = unname(mods["rna"]),
      modality_mirna        = unname(mods["mirna"]),
      modality_methylation  = unname(mods["methylation"]),
      modality_cnv          = unname(mods["cnv"]),
      modality_mutation     = unname(mods["mutation"]),
      modality_clinical     = unname(mods["clinical"]),
      n_modalities          = sum(unlist(mods[c("rna","mirna","methylation","cnv","mutation")]),
                                  na.rm = TRUE),
      has_figures           = has_figs,
      figures_dir           = if (has_figs) p$figures_dir else NA_character_,
      metrics_dir           = if (has_summ) p$metrics_dir else NA_character_,
      processed_dir         = if (has_proc) p$processed   else NA_character_,
      run_summary_path      = if (file.exists(p$run_summary)) p$run_summary else NA_character_
    )
  })
  meta <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  meta
}

build_benchmark_caches <- function(cohorts, meta) {
  # Long-format: per-cohort × model metrics
  rows <- lapply(meta[status == "completed", cohort], function(co) {
    p <- cohort_paths(co)
    s <- tryCatch(data.table::fread(p$summary_csv), error = function(e) NULL)
    if (is.null(s) || nrow(s) == 0) return(NULL)
    s[, cohort := co]
    s
  })
  long <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  setcolorder(long, c("cohort", setdiff(names(long), "cohort")))

  # Best-model table (already provided as aggregate but rebuild from long for safety)
  if (nrow(long) > 0) {
    best <- long[, .SD[which.max(c_index_test)], by = cohort,
                 .SDcols = c("model","feature_set","c_index_test","c_index_train",
                             "events_train","events_test","n_train","n_test",
                             "km_test_logrank_p")]
  } else {
    best <- data.table()
  }

  # Win counts
  if (nrow(best) > 0) {
    wins <- best[, .(wins = .N), by = model][order(-wins)]
  } else {
    wins <- data.table(model = character(), wins = integer())
  }

  saveRDS(long, file.path(CACHE_DIR, "benchmark_long.rds"))
  saveRDS(best, file.path(CACHE_DIR, "benchmark_best.rds"))
  saveRDS(wins, file.path(CACHE_DIR, "win_counts.rds"))
  list(long = long, best = best, wins = wins)
}

classify_figure <- function(fname) {
  bn <- tolower(basename(fname))
  if (grepl("^km_", bn))                    return(list(type = "km",         model = sub("^km_(.*)\\.png$","\\1", bn)))
  if (grepl("^cindex_comparison",  bn))     return(list(type = "cindex",     model = NA_character_))
  if (grepl("^cox_top_coefficients", bn))   return(list(type = "cox_coef",   model = "cox_elastic_net"))
  if (grepl("^rsf_top_features", bn))       return(list(type = "rsf_feat",   model = "random_survival_forest"))
  if (grepl("^deepsurv_training", bn))      return(list(type = "training",   model = "deepsurv"))
  if (grepl("^multibranch_training", bn))   return(list(type = "training",   model = "multibranch"))
  list(type = "other", model = NA_character_)
}

build_figure_manifest <- function(meta) {
  rows <- list()
  for (co in meta[has_figures == TRUE, cohort]) {
    fdir <- meta[cohort == co, figures_dir][1]
    if (is.na(fdir) || !dir.exists(fdir)) next
    figs <- list.files(fdir, pattern = "\\.png$", full.names = TRUE)
    for (f in figs) {
      cls <- classify_figure(f)
      rows[[length(rows) + 1L]] <- data.table(
        cohort       = co,
        figure_file  = basename(f),
        figure_path  = f,
        figure_type  = cls$type,
        model        = cls$model
      )
    }
  }
  manifest <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  saveRDS(manifest, file.path(CACHE_DIR, "figure_manifest.rds"))
  manifest
}

build_survival_cache <- function(meta) {
  out <- list()
  surv_summary_rows <- list()
  for (co in meta[status == "completed", cohort]) {
    p <- cohort_paths(co)
    y <- read_parquet_dt(file.path(p$processed, "y_survival.parquet"))
    if (is.null(y) || nrow(y) == 0) next
    setnames(y, c("OS_time","OS_event"),
             c("OS_time","OS_event"), skip_absent = TRUE)
    saveRDS(y, file.path(CACHE_DIR, "cohort_summaries", paste0(co, "_survival.rds")))
    surv_summary_rows[[co]] <- data.table(
      cohort      = co,
      n_patients  = nrow(y),
      events      = sum(y$OS_event, na.rm = TRUE),
      median_time = median(y$OS_time, na.rm = TRUE),
      max_time    = max(y$OS_time, na.rm = TRUE)
    )
    out[[co]] <- y
  }
  surv_summary <- rbindlist(surv_summary_rows, use.names = TRUE, fill = TRUE)
  saveRDS(surv_summary, file.path(CACHE_DIR, "survival_summary.rds"))
  surv_summary
}

# =============================================================================
# 2) Sample annotation + tumor/normal modality + clustering caches
# =============================================================================

# Read clinical TSV (raw GDC dump) and produce a patient-level annotation.
read_raw_clinical <- function(cohort) {
  raw_files <- cohort_raw_files(cohort)
  cf <- raw_files$clinical
  if (is.na(cf) || !file.exists(cf)) return(NULL)
  clin <- tryCatch(data.table::fread(cmd = paste("zcat", shQuote(cf))),
                   error = function(e) NULL)
  if (is.null(clin) || nrow(clin) == 0) return(NULL)

  pid_candidates <- c("submitter_id", "bcr_patient_barcode",
                      "submitter_id.samples", "_PATIENT", "sample")
  pid_col <- intersect(pid_candidates, names(clin))[1]
  if (is.na(pid_col)) return(NULL)

  age_col   <- intersect(c("age_at_initial_pathologic_diagnosis",
                           "age_at_diagnosis.diagnoses",
                           "age_at_index.demographic",
                           "age_at_diagnosis", "age_at_index"), names(clin))[1]
  sex_col   <- intersect(c("gender.demographic", "gender", "sex"), names(clin))[1]
  stage_col <- intersect(c("ajcc_pathologic_stage", "tumor_stage.diagnoses",
                           "pathologic_stage", "figo_stage"), names(clin))[1]
  vit_col   <- intersect(c("vital_status", "vital_status.demographic"), names(clin))[1]
  d2death   <- intersect(c("days_to_death"), names(clin))[1]
  d2follow  <- intersect(c("days_to_last_follow_up"), names(clin))[1]

  out <- data.table(
    patient_id = clin[[pid_col]]
  )
  out[, patient_id := tcga_patient_from_sample(patient_id)]
  if (!is.na(age_col))   out[, age          := suppressWarnings(as.numeric(clin[[age_col]]))]
  if (!is.na(sex_col))   out[, sex          := tolower(as.character(clin[[sex_col]]))]
  if (!is.na(stage_col)) out[, stage_raw    := as.character(clin[[stage_col]])]
  if (!is.na(vit_col))   out[, vital_status := tolower(as.character(clin[[vit_col]]))]
  if (!is.na(d2death))   out[, days_to_death       := suppressWarnings(as.numeric(clin[[d2death]]))]
  if (!is.na(d2follow))  out[, days_to_last_follow := suppressWarnings(as.numeric(clin[[d2follow]]))]

  # Some Xena dumps store age as negative days-from-birth
  if ("age" %in% names(out)) {
    if (sum(!is.na(out$age)) > 0 && mean(out$age < 0, na.rm = TRUE) > 0.5) {
      out[, age := round(-age / 365.25, 1)]
    }
  }
  if ("stage_raw" %in% names(out)) {
    s <- tolower(out$stage_raw)
    out[, stage := fifelse(grepl("iv", s), "IV",
                    fifelse(grepl("iii", s), "III",
                     fifelse(grepl("ii", s), "II",
                      fifelse(grepl("(\\bi\\b|stage i)", s), "I", NA_character_))))]
  }
  out <- unique(out, by = "patient_id")
  out[]
}

# Read just the column header of a (possibly gzipped) tsv to get sample IDs.
read_tsv_header <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  con <- gzfile(path, open = "r")
  on.exit(close(con), add = TRUE)
  hdr <- readLines(con, n = 1)
  if (!length(hdr)) return(character())
  strsplit(hdr, "\t", fixed = TRUE)[[1]]
}

build_sample_annotation_for_cohort <- function(cohort, meta) {
  raw <- cohort_raw_files(cohort)
  rna_path <- if (!is.null(raw$rnaseq_cnt) && !is.na(raw$rnaseq_cnt)) raw$rnaseq_cnt else raw$rnaseq_fpkm
  meth_path <- raw$methylation
  mir_path  <- raw$mirna
  cnv_path  <- raw$cnv
  mut_path  <- raw$mutation

  # Sample IDs per modality (from raw headers)
  rna_samples  <- setdiff(read_tsv_header(rna_path),  c("Ensembl_ID","sample","Hybridization REF",""))
  meth_samples <- setdiff(read_tsv_header(meth_path), c("Probe","sample",""))
  cnv_hdr <- read_tsv_header(cnv_path)
  # GDC CNV (gistic) has 3 cols per sample: <barcode>{_copy_number,_min_copy_number,_max_copy_number}
  # Pick the canonical _copy_number (NOT _min/_max) and strip the suffix.
  cn_cols <- grep("_copy_number$", cnv_hdr, value = TRUE)
  cn_cols <- cn_cols[!grepl("_(min|max)_copy_number$", cn_cols)]
  cnv_samples <- sub("_copy_number$", "", cn_cols)
  if (!length(cnv_samples)) {
    cnv_samples <- setdiff(cnv_hdr,
                           c("gene_id","gene_name","chromosome","start","end","sample",""))
  }
  # miRNA header has 3 cols per sample; we want unique sample IDs from "reads_per_million..." cols.
  mir_hdr <- read_tsv_header(mir_path)
  mir_samples <- unique(sub("^reads_per_million_miRNA_mapped_", "",
                            grep("^reads_per_million_miRNA_mapped_", mir_hdr, value = TRUE)))

  # Mutation file is long-format MAF; read the Tumor_Sample_Barcode column quickly.
  mut_samples <- character()
  if (!is.na(mut_path) && file.exists(mut_path)) {
    mut <- tryCatch(
      data.table::fread(cmd = paste("zcat", shQuote(mut_path)),
                        select = "Tumor_Sample_Barcode", showProgress = FALSE),
      error = function(e) NULL
    )
    if (!is.null(mut) && nrow(mut) > 0) mut_samples <- unique(mut$Tumor_Sample_Barcode)
  }

  # Each modality contributes (aliquot-level) sample IDs. We collapse them to
  # vial-level `sample_short` IDs (TCGA-XX-YYYY-NNV) so a single biological
  # specimen has ONE row regardless of aliquot-level duplicates across files.
  rna_short  <- tcga_sample_short(rna_samples)
  meth_short <- tcga_sample_short(meth_samples)
  cnv_short  <- tcga_sample_short(cnv_samples)
  mir_short  <- tcga_sample_short(mir_samples)
  mut_short  <- tcga_sample_short(mut_samples)

  all_short <- unique(c(rna_short, meth_short, cnv_short, mir_short, mut_short))
  all_short <- all_short[!is.na(all_short)]
  if (!length(all_short)) return(NULL)

  # Map each short id back to a representative aliquot id per modality (first hit)
  pick_alq <- function(short_target, full_short, full_full) {
    if (!length(full_short)) return(NA_character_)
    full_full[match(short_target, full_short)]
  }

  ann <- data.table(
    sample_short = all_short,
    cohort       = cohort,
    patient_id   = tcga_patient_from_sample(all_short),
    sample_type  = tcga_sample_type(all_short),
    has_rna         = all_short %in% rna_short,
    has_methylation = all_short %in% meth_short,
    has_cnv         = all_short %in% cnv_short,
    has_mirna       = all_short %in% mir_short,
    has_mutation    = all_short %in% mut_short,
    rna_aliquot     = pick_alq(all_short, rna_short,  rna_samples),
    meth_aliquot    = pick_alq(all_short, meth_short, meth_samples),
    cnv_aliquot     = pick_alq(all_short, cnv_short,  cnv_samples),
    mirna_aliquot   = pick_alq(all_short, mir_short,  mir_samples),
    mut_aliquot     = pick_alq(all_short, mut_short,  mut_samples)
  )
  ann[, sample_id := sample_short]   # canonical id used by the app

  # Attach clinical patient metadata
  clin <- read_raw_clinical(cohort)
  if (!is.null(clin)) {
    ann <- merge(ann, clin, by = "patient_id", all.x = TRUE, sort = FALSE)
  }
  # Attach survival from cohort processed y_survival (richer/cleaner)
  surv_path <- file.path(cohort_paths(cohort)$processed, "y_survival.parquet")
  if (file.exists(surv_path)) {
    y <- read_parquet_dt(surv_path)
    if (!is.null(y) && nrow(y) > 0) {
      setnames(y, "patient_id", "patient_id")
      ann <- merge(ann, y[, .(patient_id, OS_time, OS_event)],
                   by = "patient_id", all.x = TRUE, sort = FALSE)
    }
  }
  setcolorder(ann, c("sample_id","sample_short","cohort","patient_id","sample_type"))
  ann[]
}

# Collapse a (features × aliquots) numeric matrix to (features × sample_short)
# by averaging aliquot-level columns that map to the same vial-level id.
# Drops columns whose short id is NA (malformed barcodes).
collapse_aliquots_to_short <- function(m) {
  if (is.null(m) || !ncol(m)) return(m)
  short <- tcga_sample_short(colnames(m))
  keep <- !is.na(short)
  m <- m[, keep, drop = FALSE]
  short <- short[keep]
  if (!anyDuplicated(short)) {
    colnames(m) <- short
    return(m)
  }
  # Group by short id; mean across duplicate aliquots
  groups <- split(seq_along(short), short)
  out <- vapply(groups, function(idx) rowMeans(m[, idx, drop = FALSE], na.rm = TRUE),
                numeric(nrow(m)))
  rownames(out) <- rownames(m)
  out
}

build_sample_annotations <- function(cohorts_to_process, meta) {
  rows <- list()
  for (co in cohorts_to_process) {
    if (co %in% names(EXCLUDED_COHORTS)) next
    logmsg("  sample-annotation:", co)
    ann <- tryCatch(build_sample_annotation_for_cohort(co, meta),
                    error = function(e) { warning(sprintf("ann %s: %s", co, conditionMessage(e))); NULL })
    if (is.null(ann) || !nrow(ann)) next
    saveRDS(ann, file.path(CACHE_DIR, "sample_annotations", paste0(co, ".rds")))
    rows[[co]] <- ann
  }
  combined <- if (length(rows)) rbindlist(rows, use.names = TRUE, fill = TRUE) else data.table()
  saveRDS(combined, file.path(CACHE_DIR, "sample_annotations.rds"))
  combined
}

# ---------- raw modality readers (return: matrix + feature meta + sample meta)

# Generic: read a wide tsv.gz where 1st col is feature ID and remaining cols are samples.
# Returns a numeric matrix (features x samples) keeping ALL samples (tumor + normal).
read_wide_tsv_matrix <- function(path, id_col_idx = 1L, meta_cols = 0L) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  dt <- tryCatch(data.table::fread(cmd = paste("zcat", shQuote(path)), showProgress = FALSE),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  feature_ids <- dt[[id_col_idx]]
  meta <- if (meta_cols > 0L) dt[, 1:(id_col_idx + meta_cols), with = FALSE] else dt[, id_col_idx, with = FALSE]
  sample_cols <- setdiff(seq_len(ncol(dt)), seq_len(id_col_idx + meta_cols))
  if (!length(sample_cols)) return(NULL)
  m <- as.matrix(dt[, sample_cols, with = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- as.character(feature_ids)
  list(matrix = m, feature_meta = meta)
}

# Top-variance feature selection using matrixStats if available, fallback to apply.
top_variable_rows <- function(m, k) {
  if (nrow(m) <= k) return(seq_len(nrow(m)))
  v <- if (requireNamespace("matrixStats", quietly = TRUE))
    matrixStats::rowVars(m, na.rm = TRUE)
  else apply(m, 1, var, na.rm = TRUE)
  v[!is.finite(v)] <- 0
  order(v, decreasing = TRUE)[seq_len(k)]
}

tumor_normal_feature_summary <- function(m, sample_types) {
  # m: features x samples; sample_types: vector aligned to colnames(m)
  is_t <- sample_types == "tumor"
  is_n <- sample_types == "normal"
  if (sum(is_n) < 2L || sum(is_t) < 2L) return(NULL)
  if (!requireNamespace("matrixStats", quietly = TRUE)) {
    mt <- rowMeans(m[, is_t, drop = FALSE], na.rm = TRUE)
    mn <- rowMeans(m[, is_n, drop = FALSE], na.rm = TRUE)
    medt <- apply(m[, is_t, drop = FALSE], 1, median, na.rm = TRUE)
    medn <- apply(m[, is_n, drop = FALSE], 1, median, na.rm = TRUE)
    sdt <- apply(m[, is_t, drop = FALSE], 1, sd, na.rm = TRUE)
    sdn <- apply(m[, is_n, drop = FALSE], 1, sd, na.rm = TRUE)
  } else {
    mt   <- matrixStats::rowMeans2(m, cols = which(is_t), na.rm = TRUE)
    mn   <- matrixStats::rowMeans2(m, cols = which(is_n), na.rm = TRUE)
    medt <- matrixStats::rowMedians(m[, is_t, drop = FALSE], na.rm = TRUE)
    medn <- matrixStats::rowMedians(m[, is_n, drop = FALSE], na.rm = TRUE)
    sdt  <- matrixStats::rowSds(m[, is_t, drop = FALSE], na.rm = TRUE)
    sdn  <- matrixStats::rowSds(m[, is_n, drop = FALSE], na.rm = TRUE)
  }

  # Welch's two-sample t-test per row (vectorised). Cheap; we already have
  # means and SDs. Use matrixTests if available (handles NAs cleanly), else
  # roll our own welch.
  pval <- rep(NA_real_, nrow(m)); tstat <- rep(NA_real_, nrow(m))
  if (requireNamespace("matrixTests", quietly = TRUE)) {
    res <- tryCatch(
      matrixTests::row_t_welch(m[, is_t, drop = FALSE],
                                m[, is_n, drop = FALSE]),
      error = function(e) NULL
    )
    if (!is.null(res)) {
      pval  <- res$pvalue
      tstat <- res$statistic
    }
  }
  if (any(is.na(pval))) {
    nt <- sum(is_t); nn <- sum(is_n)
    se <- sqrt((sdt^2 / pmax(nt, 1L)) + (sdn^2 / pmax(nn, 1L)))
    tstat_fb <- (mt - mn) / pmax(se, 1e-12)
    df_fb <- (sdt^2/nt + sdn^2/nn)^2 /
             ((sdt^2/nt)^2/pmax(nt - 1L, 1L) + (sdn^2/nn)^2/pmax(nn - 1L, 1L))
    pval_fb <- 2 * pt(-abs(tstat_fb), df = df_fb)
    pval <- ifelse(is.na(pval), pval_fb, pval)
    tstat <- ifelse(is.na(tstat), tstat_fb, tstat)
  }
  fdr <- p.adjust(pval, method = "BH")

  out <- data.table(
    feature      = rownames(m),
    n_tumor      = sum(is_t),
    n_normal     = sum(is_n),
    mean_tumor   = mt,
    mean_normal  = mn,
    median_tumor = medt,
    median_normal= medn,
    sd_tumor     = sdt,
    sd_normal    = sdn,
    diff_mean    = mt - mn,
    log2fc_like  = mt - mn,   # for RNA/meth/mirna data is log-ish; for CNV it's a delta
    t_stat       = tstat,
    pvalue       = pval,
    fdr          = fdr,
    neglog10_p   = -log10(pmax(pval, .Machine$double.xmin))
  )
  out[order(-abs(diff_mean))]
}

precompute_pca <- function(m_top_features_x_samples, n_components = 10L) {
  # input: features x samples; we PCA the samples (so transpose)
  X <- t(m_top_features_x_samples)
  X[!is.finite(X)] <- NA
  # impute simple column means to allow prcomp
  col_means <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) {
    nas <- is.na(X[, j])
    if (any(nas)) X[nas, j] <- col_means[j]
  }
  # drop constant columns
  v <- apply(X, 2, sd, na.rm = TRUE)
  keep <- which(v > 0 & is.finite(v))
  if (length(keep) < 2L) return(NULL)
  X <- X[, keep, drop = FALSE]
  pr <- tryCatch(prcomp(X, center = TRUE, scale. = TRUE, rank. = min(n_components, ncol(X), nrow(X) - 1L)),
                 error = function(e) NULL)
  if (is.null(pr)) return(NULL)
  list(
    scores       = pr$x,
    var_explained= (pr$sdev^2) / sum(pr$sdev^2),
    n_features   = ncol(X)
  )
}

# 2D UMAP over samples on a features × samples matrix. Returns a samples × 2
# numeric matrix with rownames = sample ids, or NULL on failure.
precompute_umap <- function(m_top, n_neighbors = 15L, min_dist = 0.1) {
  if (is.null(m_top) || ncol(m_top) < 5L) return(NULL)
  if (!requireNamespace("uwot", quietly = TRUE)) return(NULL)
  X <- t(m_top); X[!is.finite(X)] <- NA
  cm <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) {
    nas <- is.na(X[, j]); if (any(nas)) X[nas, j] <- cm[j]
  }
  v <- apply(X, 2, sd, na.rm = TRUE)
  X <- X[, v > 0 & is.finite(v), drop = FALSE]
  if (ncol(X) < 2) return(NULL)
  X <- scale(X); X[!is.finite(X)] <- 0
  set.seed(42L)
  res <- tryCatch(
    uwot::umap(X,
               n_neighbors = min(n_neighbors, max(2L, nrow(X) - 1L)),
               min_dist    = min_dist,
               n_components= 2L,
               n_threads   = 1L,
               verbose     = FALSE),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)
  rownames(res) <- rownames(X)
  colnames(res) <- c("UMAP1", "UMAP2")
  res
}

write_modality_cache <- function(cohort, modality, m_top, feat_meta_top, sample_ann) {
  obj <- list(
    cohort         = cohort,
    modality       = modality,
    matrix         = m_top,             # features x samples
    feature_meta   = feat_meta_top,
    sample_meta    = sample_ann,
    created_at     = Sys.time()
  )
  saveRDS(obj, file.path(CACHE_DIR, "modalities",
                        sprintf("%s_%s_topvar.rds", cohort, modality)))
}

write_tn_summary <- function(cohort, modality, summary_dt) {
  if (is.null(summary_dt)) return(invisible(NULL))
  saveRDS(summary_dt, file.path(CACHE_DIR, "feature_summaries",
                               sprintf("%s_%s_tumor_normal.rds", cohort, modality)))
}

write_clustering_cache <- function(cohort, modality, m_top, sample_ann) {
  pca <- tryCatch(precompute_pca(m_top), error = function(e) NULL)
  if (is.null(pca)) return(invisible(NULL))
  obj <- list(
    cohort      = cohort,
    modality    = modality,
    pca_scores  = pca$scores,
    var_explained = pca$var_explained,
    sample_meta = sample_ann,
    created_at  = Sys.time()
  )
  saveRDS(obj, file.path(CACHE_DIR, "clustering",
                        sprintf("%s_%s_pca.rds", cohort, modality)))
  if (!OPT_NO_UMAP) {
    umap <- tryCatch(precompute_umap(m_top), error = function(e) NULL)
    if (!is.null(umap)) {
      uobj <- list(
        cohort      = cohort,
        modality    = modality,
        umap_scores = umap,
        sample_meta = sample_ann,
        n_neighbors = 15L,
        min_dist    = 0.1,
        created_at  = Sys.time()
      )
      saveRDS(uobj, file.path(CACHE_DIR, "clustering",
                             sprintf("%s_%s_umap.rds", cohort, modality)))
    }
  }
}

# Per-cohort modality builder (RNA / methylation / CNV / mutation / mirna).
# Reads raw files, filters to tumor+normal samples, picks top-variable features,
# writes modality cache, tumor/normal feature summary, and clustering PCA.
build_modality_for_cohort <- function(cohort, sample_ann) {
  raw <- cohort_raw_files(cohort)

  # ---------- RNA -----------------------------------------------------------
  rna_path <- if (!is.na(raw$rnaseq_cnt)) raw$rnaseq_cnt else raw$rnaseq_fpkm
  if (!is.na(rna_path) && file.exists(rna_path)) {
    logmsg("    rna   : reading", basename(rna_path))
    rna <- tryCatch(read_wide_tsv_matrix(rna_path, id_col_idx = 1L, meta_cols = 0L),
                    error = function(e) NULL)
    if (!is.null(rna)) {
      m <- rna$matrix
      mx <- suppressWarnings(max(m, na.rm = TRUE))
      if (is.finite(mx) && mx > 30) m <- log2(pmax(m, 0) + 1)
      m <- collapse_aliquots_to_short(m)
      keep_cols <- intersect(colnames(m), sample_ann$sample_id)
      if (length(keep_cols) >= 4L) {
        m <- m[, keep_cols, drop = FALSE]
        idx <- top_variable_rows(m, N_RNA)
        m_top <- m[idx, , drop = FALSE]
        ann_sub <- sample_ann[match(colnames(m_top), sample_ann$sample_id)]
        st <- ann_sub$sample_type
        feat_meta <- data.table(feature = rownames(m_top), modality = "rna")
        write_modality_cache(cohort, "rna", m_top, feat_meta, ann_sub)
        write_tn_summary(cohort, "rna", tumor_normal_feature_summary(m_top, st))
        write_clustering_cache(cohort, "rna", m_top, ann_sub)
      }
      rm(rna, m); gc(verbose = FALSE)
    }
  }

  # ---------- miRNA ---------------------------------------------------------
  if (!is.na(raw$mirna) && file.exists(raw$mirna)) {
    logmsg("    mirna : reading", basename(raw$mirna))
    mir <- tryCatch(data.table::fread(cmd = paste("zcat", shQuote(raw$mirna)),
                                      showProgress = FALSE),
                    error = function(e) NULL)
    if (!is.null(mir) && nrow(mir) > 0) {
      id_col <- names(mir)[1]
      rpm_cols <- grep("^reads_per_million_miRNA_mapped_", names(mir), value = TRUE)
      if (length(rpm_cols)) {
        sample_ids <- sub("^reads_per_million_miRNA_mapped_", "", rpm_cols)
        m <- as.matrix(mir[, rpm_cols, with = FALSE])
        storage.mode(m) <- "double"
        rownames(m) <- mir[[id_col]]
        colnames(m) <- sample_ids
        m <- collapse_aliquots_to_short(m)
        keep_cols <- intersect(colnames(m), sample_ann$sample_id)
        if (length(keep_cols) >= 4L) {
          m <- m[, keep_cols, drop = FALSE]
          m <- log2(pmax(m, 0) + 1)
          idx <- top_variable_rows(m, N_MIRNA)
          m_top <- m[idx, , drop = FALSE]
          ann_sub <- sample_ann[match(colnames(m_top), sample_ann$sample_id)]
          feat_meta <- data.table(feature = rownames(m_top), modality = "mirna")
          write_modality_cache(cohort, "mirna", m_top, feat_meta, ann_sub)
          write_tn_summary(cohort, "mirna", tumor_normal_feature_summary(m_top, ann_sub$sample_type))
          write_clustering_cache(cohort, "mirna", m_top, ann_sub)
        }
      }
      rm(mir); gc(verbose = FALSE)
    }
  }

  # ---------- CNV (gistic gene-level) --------------------------------------
  if (!is.na(raw$cnv) && file.exists(raw$cnv)) {
    logmsg("    cnv   : reading", basename(raw$cnv))
    cnv_dt <- tryCatch(data.table::fread(cmd = paste("zcat", shQuote(raw$cnv)),
                                         showProgress = FALSE),
                       error = function(e) NULL)
    if (!is.null(cnv_dt) && nrow(cnv_dt) > 0) {
      meta_cols <- intersect(c("gene_id","gene_name","chromosome","start","end","Gene Symbol"),
                             names(cnv_dt))
      # Prefer the canonical _copy_number columns when present (3-col-per-sample format).
      # `_copy_number$` also matches `_min_copy_number$` / `_max_copy_number$` so
      # we need the explicit exclusion.
      cn_cols <- grep("_copy_number$", names(cnv_dt), value = TRUE)
      cn_cols <- cn_cols[!grepl("_(min|max)_copy_number$", cn_cols)]
      if (length(cn_cols)) {
        sample_cols <- cn_cols
        new_names   <- sub("_copy_number$", "", cn_cols)
      } else {
        sample_cols <- setdiff(names(cnv_dt), meta_cols)
        new_names   <- sample_cols
      }
      if (length(sample_cols)) {
        m <- as.matrix(cnv_dt[, sample_cols, with = FALSE])
        storage.mode(m) <- "double"
        colnames(m) <- new_names
        rownames(m) <- if ("gene_name" %in% names(cnv_dt)) cnv_dt$gene_name
                       else if ("gene_id" %in% names(cnv_dt)) cnv_dt$gene_id
                       else as.character(seq_len(nrow(cnv_dt)))
        m <- collapse_aliquots_to_short(m)
        keep_cols <- intersect(colnames(m), sample_ann$sample_id)
        if (length(keep_cols) >= 4L) {
          m <- m[, keep_cols, drop = FALSE]
          idx <- top_variable_rows(m, N_CNV)
          m_top <- m[idx, , drop = FALSE]
          ann_sub <- sample_ann[match(colnames(m_top), sample_ann$sample_id)]
          fmeta_full <- if (length(meta_cols)) cnv_dt[, meta_cols, with = FALSE] else data.table(gene = rownames(m_top))
          fmeta_full <- fmeta_full[idx]
          fmeta_full[, feature := rownames(m_top)][, modality := "cnv"]
          write_modality_cache(cohort, "cnv", m_top, fmeta_full, ann_sub)
          write_tn_summary(cohort, "cnv", tumor_normal_feature_summary(m_top, ann_sub$sample_type))
          write_clustering_cache(cohort, "cnv", m_top, ann_sub)
        }
      }
      rm(cnv_dt); gc(verbose = FALSE)
    }
  }

  # ---------- Mutation (MAF -> binary patient × gene) ----------------------
  if (!is.na(raw$mutation) && file.exists(raw$mutation)) {
    logmsg("    mut   : reading", basename(raw$mutation))
    maf <- tryCatch(
      data.table::fread(cmd = paste("zcat", shQuote(raw$mutation)),
                        select = c("Hugo_Symbol","Tumor_Sample_Barcode","Variant_Classification"),
                        showProgress = FALSE),
      error = function(e) NULL
    )
    if (!is.null(maf) && nrow(maf) > 0) {
      drop_cls <- c("Silent","RNA","Intron","3'UTR","5'UTR","IGR","3'Flank","5'Flank")
      mafx <- maf[!Variant_Classification %in% drop_cls]
      gene_freq <- mafx[, .(n = uniqueN(Tumor_Sample_Barcode)), by = Hugo_Symbol][order(-n)]
      top_genes <- head(gene_freq$Hugo_Symbol, N_MUT)
      if (length(top_genes) > 0L) {
        sub <- mafx[Hugo_Symbol %in% top_genes]
        m <- dcast(sub, Hugo_Symbol ~ Tumor_Sample_Barcode,
                   value.var = "Hugo_Symbol", fun.aggregate = length)
        gene_col <- "Hugo_Symbol"
        gene_ids <- m[[gene_col]]
        m_mat <- as.matrix(m[, -1, with = FALSE])
        storage.mode(m_mat) <- "integer"
        m_mat[m_mat > 1L] <- 1L
        rownames(m_mat) <- gene_ids
        # Collapse aliquot-level mutation barcodes to vial-level (max ⇒ OR)
        short <- tcga_sample_short(colnames(m_mat))
        keep <- !is.na(short)
        m_mat <- m_mat[, keep, drop = FALSE]; short <- short[keep]
        if (anyDuplicated(short)) {
          groups <- split(seq_along(short), short)
          collapsed <- vapply(groups, function(idx) {
            as.integer(apply(m_mat[, idx, drop = FALSE], 1L, max, na.rm = TRUE))
          }, integer(nrow(m_mat)))
          rownames(collapsed) <- rownames(m_mat)
          m_mat <- collapsed
        } else {
          colnames(m_mat) <- short
        }

        # Mutation calls only exist in tumor samples - we still produce a
        # patient-level binary matrix and tag the sample type for the app.
        keep_cols <- intersect(colnames(m_mat), sample_ann$sample_id)
        if (length(keep_cols) >= 4L) {
          m_mat <- m_mat[, keep_cols, drop = FALSE]
          ann_sub <- sample_ann[match(colnames(m_mat), sample_ann$sample_id)]
          feat_meta <- gene_freq[Hugo_Symbol %in% rownames(m_mat)]
          setnames(feat_meta, "Hugo_Symbol", "feature")
          feat_meta[, modality := "mutation"]
          write_modality_cache(cohort, "mutation", m_mat, feat_meta, ann_sub)
          # tumor-vs-normal not meaningful (no normals in MAF); skip TN summary.
          # Still useful: clustering on tumor-only mutation patterns.
          tryCatch({
            obj <- list(
              cohort = cohort, modality = "mutation",
              pca_scores = NULL, var_explained = NULL,
              sample_meta = ann_sub, created_at = Sys.time(),
              note = "PCA skipped for binary mutation matrix; use Jaccard in app."
            )
            saveRDS(obj, file.path(CACHE_DIR, "clustering",
                                  sprintf("%s_%s_pca.rds", cohort, "mutation")))
          }, error = function(e) NULL)
        }
      }
      rm(maf); gc(verbose = FALSE)
    }
  }

  # ---------- Methylation (largest; gated by --no-methylation) -------------
  if (!OPT_NO_METH && !is.na(raw$methylation) && file.exists(raw$methylation)) {
    logmsg("    meth  : reading", basename(raw$methylation))
    meth <- tryCatch(read_wide_tsv_matrix(raw$methylation, id_col_idx = 1L, meta_cols = 0L),
                     error = function(e) NULL)
    if (!is.null(meth)) {
      m <- meth$matrix
      m <- collapse_aliquots_to_short(m)
      keep_cols <- intersect(colnames(m), sample_ann$sample_id)
      if (length(keep_cols) >= 4L) {
        m <- m[, keep_cols, drop = FALSE]
        # drop probes with too many NAs (>30%)
        na_frac <- rowMeans(is.na(m))
        m <- m[na_frac < 0.30, , drop = FALSE]
        idx <- top_variable_rows(m, N_METH)
        m_top <- m[idx, , drop = FALSE]
        ann_sub <- sample_ann[match(colnames(m_top), sample_ann$sample_id)]
        feat_meta <- data.table(feature = rownames(m_top), modality = "methylation")
        write_modality_cache(cohort, "methylation", m_top, feat_meta, ann_sub)
        write_tn_summary(cohort, "methylation", tumor_normal_feature_summary(m_top, ann_sub$sample_type))
        write_clustering_cache(cohort, "methylation", m_top, ann_sub)
      }
      rm(meth); gc(verbose = FALSE)
    }
  }

  invisible(TRUE)
}

# =============================================================================
# Main
# =============================================================================
main <- function() {
  logmsg("Repo:", REPO)
  logmsg("Cohort root:", COHORT_ROOT)
  logmsg("Cache dir:", CACHE_DIR)
  logmsg("Options: quick=", OPT_QUICK, " no_meth=", OPT_NO_METH,
         " no_umap=", OPT_NO_UMAP,
         " cohorts=", ifelse(is.null(OPT_COHORTS), "<all>", OPT_COHORTS))

  cohorts <- discover_cohorts()
  logmsg("Discovered cohorts:", paste(cohorts$cohort, collapse = ", "))

  # 1) Cohort metadata
  logmsg("[1/6] cohort metadata")
  meta <- build_cohort_metadata(cohorts)
  saveRDS(meta, file.path(CACHE_DIR, "cohort_metadata.rds"))
  data.table::fwrite(meta, file.path(CACHE_DIR, "cohort_metadata.csv"))

  # 2) Benchmark caches
  logmsg("[2/6] benchmark caches")
  bench <- build_benchmark_caches(cohorts, meta)

  # 3) Figure manifest
  logmsg("[3/6] figure manifest")
  figs <- build_figure_manifest(meta)

  # 4) Survival cache
  logmsg("[4/6] survival cache")
  surv <- build_survival_cache(meta)

  # 5) Sample annotations + 6) Modality caches
  if (OPT_QUICK) {
    logmsg("[5/6] skipped (--quick)")
    logmsg("[6/6] skipped (--quick)")
    sample_ann_combined <- data.table()
  } else {
    if (!is.null(OPT_COHORTS)) {
      requested <- toupper(strsplit(OPT_COHORTS, ",", fixed = TRUE)[[1]])
      requested <- ifelse(grepl("^TCGA-", requested), requested, paste0("TCGA-", requested))
      cohorts_heavy <- intersect(requested, meta[status == "completed", cohort])
    } else {
      cohorts_heavy <- meta[status == "completed", cohort]
    }
    logmsg("[5/6] sample annotations for", length(cohorts_heavy), "cohorts")
    sample_ann_combined <- build_sample_annotations(cohorts_heavy, meta)

    logmsg("[6/6] modality + tumor-normal + clustering caches")
    for (co in cohorts_heavy) {
      ann_co <- sample_ann_combined[cohort == co]
      if (!nrow(ann_co)) next
      logmsg("  cohort:", co, "samples:", nrow(ann_co),
             " tumor:", sum(ann_co$sample_type == "tumor", na.rm = TRUE),
             " normal:", sum(ann_co$sample_type == "normal", na.rm = TRUE))
      tryCatch(build_modality_for_cohort(co, ann_co),
               error = function(e) warning(sprintf("modality build failed for %s: %s",
                                                   co, conditionMessage(e))))
    }
  }

  # ------------------------------------------------------------- manifest -----
  manifest_path <- file.path(CACHE_DIR, "MANIFEST.md")
  written_files <- list.files(CACHE_DIR, recursive = TRUE)
  writeLines(c(
    "# data_cache manifest",
    "",
    paste0("Generated ", Sys.time()),
    "",
    "## Light caches (always built)",
    "- `cohort_metadata.{rds,csv}` — one row per cohort with status, sample/event sizes,",
    "  best model, modality availability and paths to figures/results.",
    "- `benchmark_long.rds` — concatenated per-cohort `summary.csv` rows for cross-cohort plots.",
    "- `benchmark_best.rds` — best model per cohort (re-derived from `benchmark_long`).",
    "- `win_counts.rds` — count of cohorts each model wins.",
    "- `figure_manifest.rds` — per-cohort PNG inventory with figure_type / model classification.",
    "- `survival_summary.rds` — per-cohort patient counts, event counts, time stats.",
    "- `cohort_summaries/<cohort>_survival.rds` — per-cohort `y_survival` table.",
    "",
    "## Sample-level caches (skipped with `--quick`)",
    "- `sample_annotations.rds` — pooled tumor/normal sample annotation across cohorts.",
    "- `sample_annotations/<cohort>.rds` — per-cohort version (recommended for app loading).",
    "",
    "## Modality caches (skipped with `--quick`)",
    "- `modalities/<cohort>_<modality>_topvar.rds` — list with top-variable feature matrix",
    "  (`features × samples`), feature metadata, sample metadata.",
    "- `feature_summaries/<cohort>_<modality>_tumor_normal.rds` — per-feature tumor-vs-normal stats",
    "  (means, medians, sds, diff_mean). Skipped for `mutation` (no normals in MAF).",
    "- `clustering/<cohort>_<modality>_pca.rds` — precomputed PCA scores + variance explained",
    "  on the same top-variable matrix. For mutation, only sample metadata is stored.",
    "- `clustering/<cohort>_<modality>_umap.rds` — precomputed 2D UMAP scores",
    "  (uwot, scaled features, n_neighbors=15, min_dist=0.1). Skipped with `--no-umap`.",
    "- Tumor-vs-normal summaries now include Welch's t-test `pvalue`, BH-corrected `fdr`",
    "  and `neglog10_p`, plus the original mean/median/sd/diff_mean columns.",
    "",
    "## Excluded cohorts",
    paste0("- `TCGA-LAML`: ", EXCLUDED_COHORTS[["TCGA-LAML"]]),
    "",
    "## Files written",
    paste0("- `", written_files, "`")
  ), manifest_path)

  # ------------------------------------------------------------- summary ------
  cat("\n========================== SUMMARY ==========================\n")
  cat(sprintf("Cohorts discovered : %d\n", nrow(cohorts)))
  cat(sprintf("  completed        : %d (%s)\n",
              nrow(meta[status == "completed"]),
              paste(meta[status == "completed", cohort], collapse = ", ")))
  cat(sprintf("  excluded         : %d (%s)\n",
              nrow(meta[status == "excluded"]),
              paste(meta[status == "excluded", cohort], collapse = ", ")))
  cat(sprintf("  partial/incomplete: %d (%s)\n",
              nrow(meta[status %in% c("partial","incomplete")]),
              paste(meta[status %in% c("partial","incomplete"), cohort], collapse = ", ")))
  cat(sprintf("Benchmark rows (long): %d\n", nrow(bench$long)))
  cat(sprintf("Figures indexed     : %d\n", nrow(figs)))
  cat(sprintf("Survival cohorts    : %d\n", nrow(surv)))
  if (exists("sample_ann_combined") && nrow(sample_ann_combined)) {
    cat(sprintf("Sample annotations  : %d samples across %d cohorts\n",
                nrow(sample_ann_combined),
                length(unique(sample_ann_combined$cohort))))
    tn <- sample_ann_combined[, .N, by = sample_type]
    cat("  by sample_type    :\n")
    for (i in seq_len(nrow(tn))) cat(sprintf("    %-8s %d\n", tn$sample_type[i], tn$N[i]))
  }
  mod_files <- list.files(file.path(CACHE_DIR, "modalities"), pattern = "\\.rds$")
  cat(sprintf("Modality caches     : %d files\n", length(mod_files)))
  tn_files  <- list.files(file.path(CACHE_DIR, "feature_summaries"), pattern = "\\.rds$")
  cat(sprintf("Tumor/normal summaries: %d files\n", length(tn_files)))
  cl_files  <- list.files(file.path(CACHE_DIR, "clustering"), pattern = "\\.rds$")
  cat(sprintf("Clustering caches   : %d files\n", length(cl_files)))
  cat(sprintf("Cache root          : %s\n", CACHE_DIR))
  cat("Manifest            :", manifest_path, "\n")
  cat("=============================================================\n")
}

main()
