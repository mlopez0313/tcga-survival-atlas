# global.R --------------------------------------------------------------------
# Loaded once at app startup. Holds:
#   - package loads (with graceful degradation),
#   - cache directory paths,
#   - small constants and lookups,
#   - eager loaders for the lightweight caches.
#
# Heavy per-cohort/per-modality data are read on demand inside the modules.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(htmltools)
  library(data.table)
  library(ggplot2)
  library(jsonlite)
})

# Optional packages -----------------------------------------------------------
.has_pkg <- function(p) suppressWarnings(requireNamespace(p, quietly = TRUE))

OPT_PKGS <- list(
  plotly         = .has_pkg("plotly"),
  DT             = .has_pkg("DT"),
  reactable      = .has_pkg("reactable"),
  shinyWidgets   = .has_pkg("shinyWidgets"),
  pheatmap       = .has_pkg("pheatmap"),
  ComplexHeatmap = .has_pkg("ComplexHeatmap"),
  circlize       = .has_pkg("circlize"),
  patchwork      = .has_pkg("patchwork"),
  viridis        = .has_pkg("viridis"),
  RColorBrewer   = .has_pkg("RColorBrewer"),
  ggrepel        = .has_pkg("ggrepel"),
  arrow          = .has_pkg("arrow"),
  fs             = .has_pkg("fs"),
  cluster        = .has_pkg("cluster"),
  factoextra     = .has_pkg("factoextra"),
  uwot           = .has_pkg("uwot"),
  survival       = .has_pkg("survival"),
  survminer      = .has_pkg("survminer"),
  matrixStats    = .has_pkg("matrixStats"),
  glue           = .has_pkg("glue"),
  dplyr          = .has_pkg("dplyr"),
  tidyr          = .has_pkg("tidyr"),
  stringr        = .has_pkg("stringr"),
  thematic       = .has_pkg("thematic")
)

if (OPT_PKGS$plotly)       suppressPackageStartupMessages(library(plotly))
if (OPT_PKGS$DT)           suppressPackageStartupMessages(library(DT))
if (OPT_PKGS$shinyWidgets) suppressPackageStartupMessages(library(shinyWidgets))
if (OPT_PKGS$dplyr)        suppressPackageStartupMessages(library(dplyr))
if (OPT_PKGS$thematic) {
  thematic::thematic_shiny(font = "auto")
}

# Paths -----------------------------------------------------------------------
# Source the config helper first so the rest of the file can rely on
# APP_PATHS / log_msg / validate_app_data.
local({
  app_dir_hint <- if (basename(getwd()) == "shiny_app") getwd() else NULL
  cfg_path <- file.path(if (!is.null(app_dir_hint)) app_dir_hint
                        else file.path(Sys.getenv("HOME"), "tcga_survival",
                                       "shiny_app"),
                        "R", "utils_config.R")
  source(cfg_path, local = FALSE)
  assign("APP_PATHS", build_app_paths(app_dir_hint), envir = globalenv())
})

# Legacy all-caps constants (kept for back-compat with existing module code)
APP_ROOT     <- APP_PATHS$app_root
CACHE_DIR    <- APP_PATHS$cache_dir
COHORT_ROOT  <- APP_PATHS$cohort_root
SUMMARY_DIR  <- APP_PATHS$summary_dir
REPO_ROOT    <- APP_PATHS$repo_root

# Constants -------------------------------------------------------------------
EXCLUDED_COHORTS <- list(
  "TCGA-LAML" = paste0(
    "TCGAbiolinks downloader could not retrieve required clinical fields ",
    "(notably `disease_response`); cohort never reached preprocessing."
  )
)

MODEL_LABELS <- c(
  cox_elastic_net        = "Cox elastic-net",
  random_survival_forest = "Random Survival Forest",
  deepsurv               = "DeepSurv",
  multibranch            = "MultiBranch (multimodal)"
)

MODALITIES <- c(
  rna         = "RNA",
  mirna       = "miRNA",
  methylation = "Methylation",
  cnv         = "CNV",
  mutation    = "Mutation"
)

DEFAULT_COHORT   <- "TCGA-LUAD"
DEFAULT_MODALITY <- "rna"

# Single source of truth for the visual language used everywhere ------------
# (utils_plots.R / module UIs read from this; CSS variables in styles.css
# are kept in sync by hand — keep the two consistent if you change colours.)
TCGA_PALETTE <- list(
  ink           = "#1f2d3d",
  ink_soft      = "#2c3e50",
  muted         = "#5a6776",
  rule          = "#dfe5ec",
  bg            = "#f7f9fb",
  card_cap      = "#eef2f6",
  accent        = "#2c7fb8",
  accent_soft   = "#9ecae1",
  warn          = "#d97706",
  good          = "#3a8b6b",
  tumor         = "#c0392b",
  normal        = "#2e86ab",
  control       = "#6c757d",
  ## qualitative palette used for clusters / models / general categories
  qual          = c("#2c7fb8","#41ab5d","#fd8d3c","#756bb1","#fa9fb5",
                    "#1d91c0","#d94801","#54278f","#bd0026","#525252"),
  ## stage gradient (I → IV)
  stage         = c(I = "#a6cee3", II = "#5d8aa8", III = "#fd8d3c", IV = "#bd0026"),
  ## model family — keeps colours stable across pages
  model_family  = c(cox_elastic_net        = "#3a8b6b",
                    random_survival_forest = "#2c7fb8",
                    deepsurv               = "#756bb1",
                    multibranch            = "#fd8d3c")
)

# Status → human label + Bootstrap badge class
status_pill <- function(status) {
  if (is.null(status) || is.na(status)) return(htmltools::tags$span(""))
  cls <- switch(as.character(status),
                completed  = "bg-success-subtle  text-success",
                excluded   = "bg-warning-subtle  text-warning",
                partial    = "bg-info-subtle     text-info",
                incomplete = "bg-secondary-subtle text-secondary",
                "bg-light text-dark")
  htmltools::tags$span(class = paste("status-pill", cls), as.character(status))
}

# Source helpers + modules ----------------------------------------------------
.r_files <- list.files(file.path(APP_ROOT, "R"), pattern = "\\.R$",
                       full.names = TRUE)
# Load utils first, then mods (alphabetical fits — utils_* sorts after mod_*,
# but R is permissive so as long as everything is loaded before server() runs
# we're fine; we explicitly source utils first for safety).
.utils <- grep("/utils_", .r_files, value = TRUE)
.mods  <- setdiff(.r_files, .utils)
for (f in .utils) source(f, local = FALSE)
for (f in .mods)  source(f, local = FALSE)

# Eager-load light caches -----------------------------------------------------
APP_DATA <- load_app_data(CACHE_DIR)

# Run startup validation. Never throws — only logs. Set
# TCGA_SURVIVAL_LOG_LEVEL=warn to silence info lines.
APP_VALIDATION <- tryCatch(
  validate_app_data(APP_PATHS, APP_DATA),
  error = function(e) {
    log_msg("validate_app_data failed: ", conditionMessage(e), level = "warn")
    NULL
  })
log_validation(APP_VALIDATION)

# Theme -----------------------------------------------------------------------
APP_THEME <- bslib::bs_theme(
  version    = 5,
  bootswatch = "minty",
  primary    = "#2c3e50",
  secondary  = "#5d8aa8",
  success    = "#3a8b6b",
  base_font  = bslib::font_google("Inter", local = FALSE),
  heading_font = bslib::font_google("Source Sans 3", local = FALSE),
  "body-bg"   = "#f7f9fb",
  "navbar-bg" = "#1f2d3d",
  "card-cap-bg" = "#eef2f6"
)
