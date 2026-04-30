# mod_overview.R --------------------------------------------------------------
# Page 1 — Overview. Editorial framing of the project, a small set of
# high-value summary cards, the cross-cohort heatmap and headline
# conclusions, then a clear "where to go next" rail.
# -----------------------------------------------------------------------------

mod_overview_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "page-header",
        div(class = "page-eyebrow", "TCGA · 18 cohorts · pan-cancer"),
        h2("A pan-cancer survival benchmark, layered with biology",
           class = "page-title"),
        p(class = "page-lede",
          "This atlas pairs a head-to-head comparison of survival models ",
          "across TCGA cohorts with ",
          tags$span(class = "lede-accent",
                    "tumor-vs-normal molecular exploration"),
          " and unsupervised structure on each modality. Use the navigation ",
          "above to follow the story from cohort landscape → biology → ",
          "structure → survival → cross-cohort patterns.")),

    layout_columns(
      col_widths = c(3,3,3,3),
      value_box(title = "Completed cohorts",
                value = textOutput(ns("kpi_cohorts"), inline = TRUE),
                showcase = bsicons::bs_icon("collection"),
                theme = "primary"),
      value_box(title = "Patients in benchmark",
                value = textOutput(ns("kpi_patients"), inline = TRUE),
                showcase = bsicons::bs_icon("people-fill"),
                theme = "secondary"),
      value_box(title = "Tumor + normal samples explored",
                value = textOutput(ns("kpi_samples"), inline = TRUE),
                showcase = bsicons::bs_icon("droplet-half"),
                theme = "success"),
      value_box(title = "Top model family",
                value = textOutput(ns("kpi_topmodel"), inline = TRUE),
                showcase = bsicons::bs_icon("trophy"),
                theme = "warning")
    ),

    layout_columns(
      col_widths = c(8, 4),
      card(card_header("Cross-cohort test C-index"),
           card_body(plotlyOutput_or_plot(ns("benchmark_heatmap"),
                                          height = "500px")),
           card_body(p(class = "text-caption",
                       "Each row is a TCGA cohort; each column a model. ",
                       "Numbers are concordance-index on the held-out test ",
                       "split. Higher = better discrimination."))),
      card(class = "right-rail",
           card_header("Headline findings"),
           card_body(uiOutput(ns("headline_text"))))
    ),

    layout_columns(
      col_widths = c(5, 7),
      card(card_header("Win count by model family"),
           card_body(plotOutput(ns("wins_plot"), height = "300px"))),
      card(card_header("Caveats & exclusions"),
           card_body(uiOutput(ns("caveats"))))
    ),

    div(class = "section-rule"),
    div(class = "page-eyebrow", style = "margin: 4px 4px 8px;",
        "Where to go next"),
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      next_step_card("ca", "Cohort atlas",
                     paste0("What data exist for each cohort, who is ",
                            "complete, and where the gaps are."),
                     icon = "grid-3x3-gap"),
      next_step_card("tn", "Tumor vs normal",
                     paste0("RNA / miRNA / methylation / CNV biology ",
                            "contrasted between tumor and matched-normal ",
                            "samples."),
                     icon = "droplet-half"),
      next_step_card("ms", "Molecular structure",
                     paste0("Find sample-level subtypes via clustering and ",
                            "ask whether they track stage or survival."),
                     icon = "diagram-3"),
      next_step_card("cc", "Cross-cohort benchmark",
                     paste0("Filter, rank, and download the cohort × model ",
                            "performance matrix."),
                     icon = "bar-chart-line")
    )
  )
}

# Helper: degrade plotlyOutput → plotOutput when plotly is missing.
plotlyOutput_or_plot <- function(id, height = "400px") {
  if (isTRUE(OPT_PKGS$plotly)) plotly::plotlyOutput(id, height = height)
  else plotOutput(id, height = height)
}

# Visual "go to page X" tile rendered on the Overview page. Clicking it
# switches the top navbar to the requested tab via `updateNavbarPage()`.
next_step_card <- function(target, title, body, icon = "arrow-right-circle") {
  tags$a(
    class = "next-step",
    href  = "#",
    onclick = sprintf(
      "Shiny.setInputValue('._goto', '%s', {priority:'event'}); return false;",
      target),
    div(class = "next-step-title",
        bsicons::bs_icon(icon),
        title),
    p(class = "next-step-body", body)
  )
}

mod_overview_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    output$kpi_cohorts <- renderText({
      m <- app_data$cohort_metadata
      if (is.null(m)) return("—")
      sprintf("%d / %d", sum(m$status == "completed"), nrow(m))
    })
    output$kpi_patients <- renderText({
      ss <- app_data$survival_summary
      if (is.null(ss)) return("—")
      formatC(sum(ss$n_patients), big.mark = ",")
    })
    output$kpi_samples <- renderText({
      sa <- read_rds_safe(file.path(app_data$cache_dir,
                                     "sample_annotations.rds"))
      if (is.null(sa) || !nrow(sa)) return("—")
      n_t <- sum(sa$sample_type == "tumor",  na.rm = TRUE)
      n_n <- sum(sa$sample_type == "normal", na.rm = TRUE)
      sprintf("%s / %s", formatC(n_t, big.mark = ","),
                          formatC(n_n, big.mark = ","))
    })
    output$kpi_topmodel <- renderText({
      w <- app_data$win_counts
      if (is.null(w) || !nrow(w)) return("—")
      sprintf("%s (%d)", model_label(w$model[1]), w$wins[1])
    })

    output$benchmark_heatmap <- if (isTRUE(OPT_PKGS$plotly))
      plotly::renderPlotly({
        p <- plot_benchmark_heatmap(app_data$benchmark_long)
        if (is.null(p)) return(plotly::plotly_empty())
        plotly::ggplotly(p, tooltip = c("x","y","fill")) |>
          plotly::config(displayModeBar = FALSE)
      })
    else
      renderPlot({ p <- plot_benchmark_heatmap(app_data$benchmark_long); print(p) })

    output$wins_plot <- renderPlot({
      p <- plot_win_counts(app_data$win_counts)
      print(p)
    })

    output$headline_text <- renderUI({
      tags$ul(class = "interp-list",
        tags$li(strong("No universal winner."),
                " Across the completed cohorts the best model family changes ",
                "by cancer type. Multimodal deep learning (MultiBranch) and ",
                "Random Survival Forest each take roughly a third of cohorts."),
        tags$li(strong("Classical baselines remain competitive."),
                " Cox elastic-net wins outright in several cohorts and is ",
                "within noise of the best model in many more."),
        tags$li(strong("DeepSurv wins selectively."),
                " It excels in cohorts with informative single-modality ",
                "signal (e.g., GBM, LGG) but rarely dominates broadly."),
        tags$li(strong("Discrimination ≠ calibration."),
                " All scores reported here are concordance index on a ",
                "single held-out test split — a ranking metric. Treat ",
                "differences below ~0.05 as inside the noise floor.")
      )
    })

    output$caveats <- renderUI({
      m <- app_data$cohort_metadata
      excl <- if (!is.null(m)) m[status == "excluded"] else data.table()
      partial <- if (!is.null(m)) m[status %in% c("partial","incomplete")] else data.table()
      tagList(
        if (nrow(excl)) tagList(
          h6(class = "caveat-h5", "Excluded cohorts"),
          tags$ul(class = "interp-list",
            lapply(seq_len(nrow(excl)), function(i) {
              tags$li(strong(excl$cohort[i]), " ", status_pill("excluded"),
                      tags$br(),
                      tags$span(class = "text-caption",
                                excl$exclusion_reason[i] %||% "no reason recorded"))
            }))) else NULL,
        if (nrow(partial)) tagList(
          h6(class = "caveat-h5", "Partial / incomplete"),
          tags$ul(class = "interp-list",
            lapply(seq_len(nrow(partial)), function(i) {
              tags$li(strong(partial$cohort[i]), " ",
                      status_pill(partial$status[i]),
                      tags$br(),
                      tags$span(class = "text-caption",
                                "Light caches still include it; heavy caches skip it."))
            }))) else NULL,
        if (!nrow(excl) && !nrow(partial)) p("No exclusions recorded.")
      )
    })
  })
}

