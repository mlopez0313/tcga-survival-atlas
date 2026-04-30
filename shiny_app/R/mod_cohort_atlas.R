# mod_cohort_atlas.R ----------------------------------------------------------
# Page 2 — Cohort Atlas. What data exist, modality availability, sizes.
# 3-pane layout with a contextual right-rail per selected cohort.
# -----------------------------------------------------------------------------

mod_cohort_atlas_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "page-header",
        div(class = "page-eyebrow", "Inventory · 18 TCGA cohorts"),
        h2("Cohort atlas", class = "page-title"),
        p(class = "page-lede",
          "What data are available for each TCGA cohort, which cohorts ",
          "made it through the modeling pipeline, and where the gaps are. ",
          "Click any row to see the cohort detail in the right rail.")),
    layout_sidebar(
      sidebar = sidebar(
        width = 280, open = "always",
        h6("Filters"),
        checkboxInput(ns("only_completed"), "Show completed cohorts only", FALSE),
        sliderInput(ns("min_n"), "Minimum total samples",
                    min = 0, max = 800, value = 0, step = 10),
        checkboxGroupInput(ns("modalities_required"),
                           "Must include modality",
                           choices = MODALITIES,
                           selected = character(0)),
        hr(),
        downloadButton(ns("dl_csv"), "Download cohort metadata (CSV)",
                       class = "btn-sm btn-outline-secondary")
      ),

      layout_columns(
        col_widths = c(8, 4),

        # ----------------------------------------------------- center -----
        tagList(
          layout_columns(
            col_widths = c(7, 5),
            card(card_header("Modality availability"),
                 card_body(plotOutput(ns("avail_heatmap"), height = "520px"))),
            card(card_header("Train / test sizes"),
                 card_body(plotOutput(ns("sample_bars"), height = "240px")),
                 card_body(plotOutput(ns("tn_bars"), height = "240px")))
          ),
          card(card_header("Cohort metadata"),
               card_body(DT_or_table_ui(ns("meta_table"))))
        ),

        # ----------------------------------------------------- right rail -
        tagList(
          card(class = "right-rail",
               card_header("Cohort details"),
               card_body(uiOutput(ns("details")))),
          card(class = "right-rail",
               card_header("Status legend"),
               card_body(
                 tags$ul(class = "interp-list",
                   tags$li(status_pill("completed"),
                           " preprocessing + modeling finished, full caches available."),
                   tags$li(status_pill("partial"),
                           " preprocessing finished but some artefacts missing."),
                   tags$li(status_pill("incomplete"),
                           " download done but preprocessing did not complete; ",
                           "in light caches only."),
                   tags$li(status_pill("excluded"),
                           " removed from the benchmark; reason recorded in metadata."))))
        )
      )
    )
  )
}

mod_cohort_atlas_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {

    selected_cohort <- reactiveVal(DEFAULT_COHORT)

    filtered <- reactive({
      m <- app_data$cohort_metadata
      if (is.null(m)) return(NULL)
      d <- as.data.table(copy(m))
      if (isTRUE(input$only_completed)) d <- d[status == "completed"]
      if (length(input$min_n) && input$min_n > 0)
        d <- d[(coalesce_int(n_train) + coalesce_int(n_test)) >= input$min_n]
      if (length(input$modalities_required)) {
        for (mod in input$modalities_required) {
          col <- paste0("modality_", mod)
          if (col %in% names(d)) d <- d[d[[col]] == TRUE]
        }
      }
      d
    })

    output$avail_heatmap <- renderPlot({
      print(plot_modality_heatmap(filtered()))
    })

    output$sample_bars <- renderPlot({
      d <- filtered(); if (is.null(d) || !nrow(d)) return(invisible(NULL))
      d <- d[status == "completed"]
      if (!nrow(d)) return(invisible(NULL))
      df <- data.table(
        cohort = factor(d$cohort, levels = d$cohort[order(-(coalesce_int(d$n_train)+coalesce_int(d$n_test)))]),
        train  = coalesce_int(d$n_train),
        test   = coalesce_int(d$n_test)
      )
      long <- melt(df, id.vars = "cohort", variable.name = "split",
                   value.name = "n")
      ggplot2::ggplot(long, ggplot2::aes(x = cohort, y = n, fill = split)) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::scale_fill_manual(values = c(train = TCGA_PALETTE$accent,
                                              test  = TCGA_PALETTE$tumor),
                                    name = NULL) +
        ggplot2::labs(x = NULL, y = "n patients",
                      title = "Train / test split sizes") +
        theme_editorial() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
    })

    output$tn_bars <- renderPlot({
      cohorts <- if (!is.null(filtered())) filtered()$cohort else character()
      rows <- lapply(cohorts, function(co) {
        ann <- load_sample_annotation(co)
        if (is.null(ann) || !nrow(ann)) return(NULL)
        data.table(cohort = co,
                   tumor  = sum(ann$sample_type == "tumor",  na.rm = TRUE),
                   normal = sum(ann$sample_type == "normal", na.rm = TRUE))
      })
      df <- rbindlist(rows, fill = TRUE)
      if (!nrow(df)) return(invisible(NULL))
      df[, cohort := factor(cohort, levels = df$cohort[order(-(tumor + normal))])]
      long <- melt(df, id.vars = "cohort", variable.name = "type",
                   value.name = "n")
      long[, type := factor(type, levels = c("tumor","normal"))]
      ggplot2::ggplot(long, ggplot2::aes(x = cohort, y = n, fill = type)) +
        ggplot2::geom_col(width = 0.7) +
        tn_fill_scale(name = NULL) +
        ggplot2::labs(x = NULL, y = "n samples",
                      title = "Tumor vs normal samples (vial-level barcodes)") +
        theme_editorial() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
    })

    if (isTRUE(OPT_PKGS$DT)) {
      output$meta_table <- DT::renderDT({
        d <- filtered(); if (is.null(d)) return(NULL)
        cols <- intersect(c("cohort","status","n_train","n_test",
                            "events_train","events_test","best_model",
                            "best_test_cindex","n_modalities"),
                          names(d))
        out <- d[, ..cols]
        if ("best_model" %in% names(out))
          out[, best_model := model_label(best_model)]
        dt <- DT::datatable(out, rownames = FALSE,
                      selection = list(mode = "single", selected = 1),
                      options = list(pageLength = 12, dom = "ftip",
                                     scrollX = TRUE),
                      class = "compact stripe")
        if ("best_test_cindex" %in% names(out))
          dt <- DT::formatRound(dt, "best_test_cindex", digits = 3)
        dt
      })
      observeEvent(input$meta_table_rows_selected, {
        d <- filtered(); if (is.null(d) || !nrow(d)) return()
        sel <- input$meta_table_rows_selected
        if (length(sel)) selected_cohort(as.character(d$cohort[sel[1]]))
      })
    } else {
      output$meta_table <- renderTable({
        d <- filtered(); if (is.null(d)) return(NULL); d
      })
    }

    output$details <- renderUI({
      co <- selected_cohort()
      m  <- app_data$cohort_metadata
      if (is.null(m)) return(empty_state("No metadata cache",
                                         "Run scripts/prepare_app_data.R first.",
                                         icon = "exclamation-circle"))
      row <- m[cohort == co]
      if (!nrow(row)) return(empty_state("Select a cohort row",
                                         "Click any row in the table to see its details here."))
      ann <- load_sample_annotation(co)
      n_t <- if (!is.null(ann)) sum(ann$sample_type == "tumor", na.rm = TRUE) else NA
      n_n <- if (!is.null(ann)) sum(ann$sample_type == "normal", na.rm = TRUE) else NA
      mod_present <- names(MODALITIES)[
        sapply(names(MODALITIES), function(mn) {
          col <- paste0("modality_", mn); isTRUE(row[[col]])
        })]
      tagList(
        div(style = "display:flex; align-items:center; gap:10px;",
            h4(co, style = "margin:0;"),
            status_pill(row$status)),
        if (!is.na(row$exclusion_reason))
          tags$blockquote(row$exclusion_reason)
        else NULL,
        tags$dl(class = "summary-dl",
          tags$dt("Patients (train/test)"),
            tags$dd(sprintf("%s / %s",
                            row$n_train %||% "—", row$n_test %||% "—")),
          tags$dt("Events (train/test)"),
            tags$dd(sprintf("%s / %s",
                            row$events_train %||% "—",
                            row$events_test %||% "—")),
          tags$dt("Best model"),
            tags$dd(if (!is.na(row$best_model))
                      sprintf("%s · C = %.3f",
                              model_label(row$best_model),
                              row$best_test_cindex)
                    else "—"),
          tags$dt("Tumor / normal"),
            tags$dd(sprintf("%s / %s",
                            if (is.na(n_t)) "—" else n_t,
                            if (is.na(n_n)) "—" else n_n)),
          tags$dt("Modalities"),
            tags$dd(if (length(mod_present))
                      paste(modality_label(mod_present), collapse = ", ")
                    else "—")
        )
      )
    })

    output$dl_csv <- downloadHandler(
      filename = function() "cohort_metadata.csv",
      content  = function(file) {
        d <- filtered(); if (is.null(d)) d <- app_data$cohort_metadata
        data.table::fwrite(d, file)
      }
    )
  })
}

coalesce_int <- function(x) {
  v <- suppressWarnings(as.integer(x))
  v[is.na(v)] <- 0L; v
}
