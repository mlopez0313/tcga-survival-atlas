# mod_cross_cohort.R ----------------------------------------------------------
# Page 6 — Cross-cohort benchmark deep dive.
# 3-pane layout: filters | heatmap+wins+scatter+table | interpretation rail.
# -----------------------------------------------------------------------------

mod_cross_cohort_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "page-header",
        div(class = "page-eyebrow", "Pan-cancer · 18 cohorts × 4 models"),
        h2("Cross-cohort benchmark", class = "page-title"),
        p(class = "page-lede",
          "Filter by event count, sample size or model family, and watch ",
          "the heatmap, win counts and performance-vs-events plot update ",
          "in lockstep.")),
    layout_sidebar(
      sidebar = sidebar(
        width = 280, open = "always",
        h6("Filters"),
        sliderInput(ns("min_events"), "Min total events",
                    min = 0, max = 400, value = 0, step = 5),
        sliderInput(ns("min_n"), "Min total samples",
                    min = 0, max = 800, value = 0, step = 10),
        checkboxGroupInput(ns("models"), "Model families",
                           choices = setNames(names(MODEL_LABELS), MODEL_LABELS),
                           selected = names(MODEL_LABELS)),
        hr(),
        downloadButton(ns("dl_csv"), "Download benchmark (CSV)",
                       class = "btn-sm btn-outline-secondary")
      ),
      layout_columns(
        col_widths = c(8, 4),

        # ----------------------------------------------------- center -----
        tagList(
          card(card_header("Cohort × model heatmap"),
               card_body(plotlyOutput_or_plot(ns("heatmap"), height = "560px")),
               card_body(p(class = "text-caption",
                           "Hover any cell for the underlying numbers; the ",
                           "color scale is centered at C = 0.65 (typical ",
                           "discrimination threshold)."))),
          layout_columns(
            col_widths = c(5, 7),
            card(card_header("Win count by model"),
                 card_body(plotOutput(ns("wins"), height = "260px"))),
            card(card_header("Performance vs events"),
                 card_body(plotOutput(ns("perf_vs_size"), height = "260px")))
          ),
          card(card_header("Sortable benchmark table"),
               card_body(DT_or_table_ui(ns("table"))))
        ),

        # ----------------------------------------------------- right rail -
        tagList(
          card(class = "right-rail",
               card_header("Interpretation"),
               card_body(uiOutput(ns("interp")))),
          card(class = "right-rail",
               card_header("Reading the heatmap"),
               card_body(
                 p("Each row is a cohort, each column a model. The number ",
                   "in each cell is the test C-index on the held-out split."),
                 p("Look for: (a) cohorts where one model clearly wins, ",
                   "(b) cohorts where all four cluster within ±0.05 — the ",
                   "noise floor for that cohort."),
                 p(class = "text-caption",
                   "C ≈ 0.5 is uninformative; C ≥ 0.7 is generally usable; ",
                   "C ≥ 0.8 is strong but rare in pan-cancer benchmarks.")))
        )
      )
    )
  )
}

DT_or_table_ui <- function(id) {
  if (isTRUE(OPT_PKGS$DT)) DT::DTOutput(id) else tableOutput(id)
}

mod_cross_cohort_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {

    filtered <- reactive({
      bl <- app_data$benchmark_long
      if (is.null(bl)) return(NULL)
      d <- as.data.table(bl)
      if (length(input$models)) d <- d[model %in% input$models]
      d <- d[(events_train + events_test) >= input$min_events]
      d <- d[(n_train + n_test) >= input$min_n]
      d
    })

    output$heatmap <- if (isTRUE(OPT_PKGS$plotly))
      plotly::renderPlotly({
        d <- filtered(); p <- plot_benchmark_heatmap(d)
        if (is.null(p)) return(plotly::plotly_empty())
        plotly::ggplotly(p, tooltip = c("x","y","fill")) |>
          plotly::config(displayModeBar = FALSE)
      })
    else renderPlot({ p <- plot_benchmark_heatmap(filtered()); print(p) })

    output$wins <- renderPlot({
      d <- filtered()
      if (is.null(d) || !nrow(d)) return(invisible(NULL))
      best <- d[, .SD[which.max(c_index_test)], by = cohort]
      wc <- best[, .(wins = .N), by = model][order(-wins)]
      print(plot_win_counts(wc))
    })

    output$perf_vs_size <- renderPlot({
      d <- filtered(); p <- plot_perf_vs_size(d); print(p)
    })

    if (isTRUE(OPT_PKGS$DT)) {
      output$table <- DT::renderDT({
        d <- filtered(); if (is.null(d)) return(NULL)
        cols <- intersect(c("cohort","model","feature_set","c_index_train",
                            "c_index_test","n_train","n_test",
                            "events_train","events_test"), names(d))
        out <- d[, ..cols][order(cohort, -c_index_test)]
        out[, model := model_label(model)]
        dt <- DT::datatable(out, rownames = FALSE,
                      options = list(pageLength = 12, dom = "ftip",
                                     scrollX = TRUE),
                      class = "compact stripe")
        rd <- intersect(c("c_index_train","c_index_test"), names(out))
        if (length(rd)) dt <- DT::formatRound(dt, rd, digits = 3)
        dt
      })
    } else {
      output$table <- renderTable({
        d <- filtered(); if (is.null(d)) return(NULL)
        head(d[order(cohort, -c_index_test)], 30L)
      })
    }

    output$interp <- renderUI({
      d <- filtered()
      if (is.null(d) || !nrow(d))
        return(empty_state("No cohorts match", "Loosen the filters to see results."))
      best <- d[, .SD[which.max(c_index_test)], by = cohort]
      wc   <- best[, .(wins = .N), by = model][order(-wins)]
      n_co <- length(unique(d$cohort))
      tight <- d[, .(spread = max(c_index_test) - min(c_index_test)),
                 by = cohort][spread < 0.05]
      tags$ul(class = "interp-list",
        tags$li(strong(sprintf("%d cohort%s in view.", n_co,
                                if (n_co == 1) "" else "s"))),
        if (nrow(wc))
          tags$li(strong("Wins under filters: "),
                  paste(sprintf("%s (%d)", model_label(wc$model), wc$wins),
                        collapse = " · "))
        else NULL,
        if (nrow(tight))
          tags$li(strong("Noise-floor cohorts (Δ C < 0.05): "),
                  paste(tight$cohort, collapse = ", "), ".")
        else NULL,
        tags$li(strong("No universal best model is observed."),
                " The best model family depends on the cancer type."),
        tags$li(strong("Cox elastic-net stays competitive"),
                " — a strong baseline before reaching for deep models."),
        tags$li(strong("DeepSurv wins selectively"),
                " — especially in cohorts with informative single-modality ",
                "signal."),
        tags$li(strong("Trust the small print."),
                " For cohorts with < 30 test events, treat differences ",
                "below 0.05 in C-index as inside the noise floor.")
      )
    })

    output$dl_csv <- downloadHandler(
      filename = function() "tcga_survival_benchmark.csv",
      content  = function(file) {
        d <- filtered(); if (is.null(d)) d <- app_data$benchmark_long
        data.table::fwrite(d, file)
      }
    )
  })
}
