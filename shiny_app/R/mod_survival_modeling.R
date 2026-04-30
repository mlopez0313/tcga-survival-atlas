# mod_survival_modeling.R -----------------------------------------------------
# Page 5 — Per-cohort survival modeling: metrics + saved figures.
# 3-pane layout: cohort/model picker | metrics + saved figures | interpretation.
# -----------------------------------------------------------------------------

mod_survival_modeling_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "page-header",
        div(class = "page-eyebrow", "Per-cohort modeling · saved figures"),
        h2("Survival modeling", class = "page-title"),
        p(class = "page-lede",
          "Inspect within-cohort survival benchmarks: train/test C-index, ",
          "Kaplan–Meier curves, training curves and feature importance — ",
          "all rendered from the saved figures the modeling pipeline ",
          "wrote to disk.")),
    layout_sidebar(
      sidebar = sidebar(
        width = 280, open = "always",
        h6("Cohort"),
        selectInput(ns("cohort"), NULL, choices = NULL, selected = DEFAULT_COHORT),
        hr(),
        h6("Models"),
        checkboxGroupInput(ns("models"), NULL,
                           choices = setNames(names(MODEL_LABELS), MODEL_LABELS),
                           selected = names(MODEL_LABELS)),
        hr(),
        downloadButton(ns("dl_summary"), "Download cohort summary (CSV)",
                       class = "btn-sm btn-outline-secondary")
      ),
      layout_columns(
        col_widths = c(8, 4),

        # ----------------------------------------------------- center -----
        tagList(
          layout_columns(
            col_widths = c(7, 5),
            card(card_header("Per-model metric summary"),
                 card_body(DT_or_table_ui(ns("summary_table")))),
            card(card_header("Train vs test C-index"),
                 card_body(plotOutput(ns("ttest"), height = "320px")))
          ),
          card(card_header("Saved figures"),
               card_body(uiOutput(ns("figures"))))
        ),

        # ----------------------------------------------------- right rail -
        tagList(
          card(class = "right-rail",
               card_header("Interpretation"),
               card_body(uiOutput(ns("interp")))),
          card(class = "right-rail",
               card_header("Reading the panels"),
               card_body(
                 p("Train vs test bars: a large train ↔ test gap signals ",
                   "overfitting; a small gap with low test C signals ",
                   "structural underperformance."),
                 p("KM curves: visual separation only — formal evidence ",
                   "comes from the log-rank p in each figure caption."),
                 p(class = "text-caption",
                   "Cohorts with < 30 test events have very wide CIs; treat ",
                   "model-to-model differences as a soft preference here.")))
        )
      )
    )
  )
}

mod_survival_modeling_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {

    observe({
      m <- app_data$cohort_metadata
      choices <- if (!is.null(m))
        as.character(m[status == "completed", cohort]) else character()
      updateSelectInput(session, "cohort", choices = choices,
                        selected = isolate(input$cohort) %||% DEFAULT_COHORT)
    })

    cohort_summary <- reactive({
      co <- input$cohort; req(co)
      d <- load_cohort_summary_csv(co)
      if (is.null(d)) {
        bl <- app_data$benchmark_long
        if (!is.null(bl)) d <- as.data.table(bl)[cohort == co]
      }
      if (is.null(d) || !nrow(d)) return(NULL)
      if (length(input$models)) d <- d[model %in% input$models]
      d
    })

    if (isTRUE(OPT_PKGS$DT)) {
      output$summary_table <- DT::renderDT({
        d <- cohort_summary(); if (is.null(d)) return(NULL)
        cols <- intersect(c("model","feature_set","c_index_train",
                            "c_index_test","km_test_logrank_p",
                            "n_train","n_test","events_train","events_test"),
                          names(d))
        out <- d[, ..cols]
        if ("model" %in% names(out)) out[, model := model_label(model)]
        dt <- DT::datatable(out, rownames = FALSE,
                      options = list(pageLength = 10, dom = "tip",
                                     scrollX = TRUE),
                      class = "compact stripe")
        round_cols <- intersect(c("c_index_train","c_index_test"), names(out))
        if (length(round_cols)) dt <- DT::formatRound(dt, round_cols, digits = 3)
        sig_cols <- intersect("km_test_logrank_p", names(out))
        if (length(sig_cols)) dt <- DT::formatSignif(dt, sig_cols, digits = 3)
        dt
      })
    } else {
      output$summary_table <- renderTable({
        d <- cohort_summary(); if (is.null(d)) return(NULL); d
      })
    }

    output$ttest <- renderPlot({
      d <- cohort_summary()
      if (is.null(d) || !nrow(d)) return(invisible(NULL))
      df <- data.table(model = model_label(d$model),
                       train = d$c_index_train,
                       test  = d$c_index_test)
      long <- melt(df, id.vars = "model", variable.name = "split",
                   value.name = "c_index")
      ggplot2::ggplot(long, ggplot2::aes(x = model, y = c_index,
                                          fill = split)) +
        ggplot2::geom_col(position = "dodge", width = 0.65) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", c_index)),
                           position = ggplot2::position_dodge(width = 0.65),
                           vjust = -0.4, size = 3,
                           color = TCGA_PALETTE$ink_soft) +
        ggplot2::scale_fill_manual(values = c(train = TCGA_PALETTE$accent,
                                              test  = TCGA_PALETTE$tumor),
                                    name = NULL) +
        ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                            color = "#b0b8c1") +
        ggplot2::scale_y_continuous(limits = c(0.4, 1),
                                     breaks = seq(0.4, 1, 0.1)) +
        ggplot2::labs(x = NULL, y = "C-index",
                      title = sprintf("Train vs test — %s", input$cohort),
                      subtitle = "0.5 dashed line = no discrimination") +
        theme_editorial() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
    })

    output$figures <- renderUI({
      co <- input$cohort; req(co)
      figs <- list_cohort_figures(co)
      if (!length(figs))
        return(empty_state("No saved figures",
                           sprintf("The modeling pipeline did not write any figures for %s.",
                                   co),
                           icon = "image"))
      tagList(
        layout_columns(col_widths = rep_len(6, length(figs)),
          !!! lapply(figs, function(f) {
            tryCatch({
              raw <- base64enc::base64encode(f)
              card(card_header(basename(f)),
                   card_body(tags$img(
                     src = paste0("data:image/png;base64,", raw),
                     style = "max-width: 100%; height: auto;")))
            }, error = function(e) card(card_header(basename(f)),
                                        card_body(empty_state(
                                          "Could not load image",
                                          conditionMessage(e),
                                          icon = "x-circle"))))
          })
        )
      )
    })

    output$interp <- renderUI({
      d <- cohort_summary()
      if (is.null(d) || !nrow(d))
        return(empty_state("No metric summary",
                           "Pick a different cohort or run the modeling pipeline first.",
                           icon = "info-circle"))
      msgs <- interp_survival_cohort(d, cohort = input$cohort)
      tags$ul(class = "interp-list",
              lapply(msgs, function(m) tags$li(m$children)))
    })

    output$dl_summary <- downloadHandler(
      filename = function() sprintf("%s_summary.csv", input$cohort),
      content  = function(file) {
        d <- cohort_summary(); if (is.null(d)) d <- data.table()
        data.table::fwrite(d, file)
      }
    )
  })
}
