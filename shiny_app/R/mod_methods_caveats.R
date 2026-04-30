# mod_methods_caveats.R -------------------------------------------------------
# Page 7 — Methods & Caveats. Pulls in repo README/METHODS/RESULTS_SUMMARY.
# -----------------------------------------------------------------------------

mod_methods_caveats_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "page-header",
        div(class = "page-eyebrow", "Methodology · caveats · provenance"),
        h2("Methods & caveats", class = "page-title"),
        p(class = "page-lede",
          "How the data was prepared, which models were trained, the ",
          "exclusions that matter, and the standard interpretive caveats ",
          "for any survival benchmark.")),
    layout_columns(
      col_widths = c(7, 5),
      card(card_header("Methodology summary"),
           card_body(uiOutput(ns("methods")))),
      card(class = "right-rail",
           card_header("Caveats"),
           card_body(uiOutput(ns("caveats"))))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(card_header("Repo files"),
           card_body(uiOutput(ns("repo_files")))),
      card(card_header("Cache manifest"),
           card_body(uiOutput(ns("manifest"))))
    ),
    card(card_header(tagList(
           bsicons::bs_icon("activity"), " Diagnostics & runtime checks")),
         card_body(uiOutput(ns("diag_summary"))),
         card_body(DT_or_table_ui(ns("diag_table"))),
         card_body(p(class = "text-caption",
                     "Run from the shell with ",
                     tags$code("Rscript shiny_app/scripts/sanity_check.R"),
                     " for a non-interactive equivalent. Use environment ",
                     "variables (",
                     tags$code("TCGA_SURVIVAL_CACHE_DIR"), ", ",
                     tags$code("TCGA_SURVIVAL_COHORT_ROOT"),
                     ") to relocate paths for deployment.")))
  )
}

mod_methods_caveats_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    output$methods <- renderUI({
      tagList(
        h5("Pipeline"),
        p("Per-cohort GDC data are pulled with TCGAbiolinks / Xena ",
          "harmonized TSVs, then preprocessed into patient-level matrices: ",
          "clinical, RNA (htseq counts), miRNA, methylation 450K, GISTIC ",
          "CNV and a MAF-derived binary mutation matrix. Modeling pipelines ",
          "train Cox elastic-net, Random Survival Forest, DeepSurv, and a ",
          "multimodal MultiBranch network. Performance is reported as test ",
          "concordance (C-index) on a single held-out split per cohort."),
        h5("Tumor vs normal exploration"),
        p("The modeling pipeline filters to primary tumors. For tumor-vs-",
          "normal exploration the app reads the ", em("raw"),
          " molecular files and re-derives sample-level matrices, keeping ",
          "all sample types. Sample IDs are normalized to vial-level ",
          tags$code("TCGA-XX-YYYY-NNV"),
          " barcodes; aliquot-level numeric duplicates are mean-collapsed, ",
          "mutation duplicates OR-collapsed."),
        h5("Statistical tests"),
        p("Tumor vs normal: Welch's two-sample t-test per feature on the ",
          "top-variable subset (1000 RNA / 200 miRNA / 2000 methylation ",
          "probes / 1000 CNV genes / 300 mutated genes by default), with ",
          "BH-corrected FDR alongside raw p-values. Cluster vs survival: ",
          "log-rank test on Kaplan–Meier curves. Cluster vs stage: χ² ",
          "test on the cluster × stage contingency table."),
        h5("Clustering"),
        p("Per-cohort/modality clustering uses the same top-variable ",
          "feature subset as above. PCA scores and UMAP coordinates ",
          "(scaled features, n_neighbors = 15, min_dist = 0.1) are ",
          "precomputed; k-means and hierarchical (Ward.D2) are exposed in ",
          "the app and re-run on every parameter change.")
      )
    })

    output$caveats <- renderUI({
      tags$ul(class = "interp-list",
        tags$li(strong("LAML omission. "),
                "TCGAbiolinks could not pull required clinical fields ",
                "(notably ", tags$code("disease_response"),
                ") for ", strong("TCGA-LAML"),
                "; that cohort is excluded from the benchmark and shown ",
                "explicitly across the app rather than dropped silently."),
        tags$li(strong("Discrimination ≠ calibration. "),
                "All scores are concordance index — a ranking metric. ",
                "Predicted risk magnitudes are not calibrated and should ",
                "not be interpreted as absolute hazards."),
        tags$li(strong("Small cohorts / few events. "),
                "CHOL, DLBC, KICH have < 10 test events; treat ranking ",
                "differences below ~0.05 in C-index as inside the noise floor."),
        tags$li(strong("Tumor-only models, tumor + normal exploration. "),
                "All survival models were trained on tumor samples only. ",
                "The tumor-vs-normal views are descriptive and do not feed ",
                "back into model training."),
        tags$li(strong("Cluster labels are unsupervised. "),
                "They reflect the dominant variance in the chosen modality ",
                "and feature subset. Treat them as hypotheses to validate ",
                "with stage / molecular subtype rather than ground truth."),
        tags$li(strong("Multiple-testing. "),
                "The volcano FDR controls the false-discovery rate within ",
                "the displayed top-variable subset; it does not generalize ",
                "to the full transcriptome / methylome.")
      )
    })

    output$repo_files <- renderUI({
      paths <- c(
        README           = file.path(APP_PATHS$repo_root,    "README.md"),
        METHODS          = file.path(APP_PATHS$repo_root,    "METHODS.md"),
        RESULTS_SUMMARY  = file.path(APP_PATHS$repo_root,    "RESULTS_SUMMARY.md"),
        AGGREGATE_README = file.path(APP_PATHS$summary_dir, "README.md")
      )
      out <- list()
      for (nm in names(paths)) {
        body <- read_text_safe(paths[[nm]])
        out[[length(out) + 1L]] <- tags$details(
          tags$summary(strong(nm), tags$code(paths[[nm]])),
          if (nzchar(body)) tags$pre(body)
          else empty_state("File not found",
                            sprintf("No file at %s.", paths[[nm]]),
                            icon = "file-earmark"))
      }
      tagList(out)
    })

    output$manifest <- renderUI({
      mp <- app_data$manifest_md
      if (!file.exists(mp))
        return(empty_state("No cache manifest yet",
                           "Run scripts/prepare_app_data.R to generate one.",
                           icon = "file-earmark-text"))
      tags$pre(read_text_safe(mp))
    })

    # ---- Diagnostics ------------------------------------------------
    diag_results <- reactive({
      tryCatch(validate_app_data(APP_PATHS, APP_DATA),
               error = function(e) {
                 data.table::data.table(
                   check = "validation",
                   status = "missing",
                   detail = conditionMessage(e))
               })
    })

    output$diag_summary <- renderUI({
      v <- diag_results()
      if (is.null(v) || !nrow(v))
        return(empty_state("Validator returned nothing",
                           "Re-run with TCGA_SURVIVAL_LOG_LEVEL=debug to investigate.",
                           icon = "exclamation-circle"))
      n_ok      <- sum(v$status == "ok",      na.rm = TRUE)
      n_warn    <- sum(v$status == "warn",    na.rm = TRUE)
      n_missing <- sum(v$status == "missing", na.rm = TRUE)
      tagList(
        tags$dl(class = "summary-dl",
          tags$dt("App root"),    tags$dd(tags$code(APP_PATHS$app_root)),
          tags$dt("Cache dir"),   tags$dd(tags$code(APP_PATHS$cache_dir)),
          tags$dt("Cohort root"), tags$dd(tags$code(APP_PATHS$cohort_root)),
          tags$dt("Summary dir"), tags$dd(tags$code(APP_PATHS$summary_dir)),
          tags$dt("Status"),
            tags$dd(status_pill(if (n_missing) "incomplete"
                                else if (n_warn) "partial"
                                else "completed"),
                    sprintf(" %d ok · %d warn · %d missing",
                            n_ok, n_warn, n_missing))
        ),
        if (n_missing)
          tags$blockquote(
            "Missing items below typically mean the data prep step has not ",
            "been run yet, or the cohort root is not mounted. Run ",
            tags$code("Rscript shiny_app/scripts/prepare_app_data.R"),
            " from the repo root.")
        else NULL
      )
    })

    if (isTRUE(OPT_PKGS$DT)) {
      output$diag_table <- DT::renderDT({
        v <- diag_results(); if (is.null(v)) return(NULL)
        DT::datatable(v, rownames = FALSE,
                      options = list(pageLength = 15, dom = "ftip",
                                     scrollX = TRUE),
                      class = "compact stripe") |>
          DT::formatStyle("status",
                           target = "cell",
                           backgroundColor = DT::styleEqual(
                             c("ok","warn","missing"),
                             c("#dff5e6","#fff3cd","#fde2e1")),
                           fontWeight = "600")
      })
    } else {
      output$diag_table <- renderTable({
        v <- diag_results(); if (is.null(v)) return(NULL); v
      })
    }
  })
}
