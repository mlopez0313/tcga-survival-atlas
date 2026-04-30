# mod_molecular_structure.R ---------------------------------------------------
# Page 4 — Molecular Structure: clustering + dim-red + clinical/survival
# interrogation by cluster. 3-pane layout with a contextual right rail.
# -----------------------------------------------------------------------------

mod_molecular_structure_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "page-header",
        div(class = "page-eyebrow", "Unsupervised structure · per cohort"),
        h2("Molecular structure", class = "page-title"),
        p(class = "page-lede",
          "Find sample-level structure within a cohort × modality and ",
          "interrogate the resulting clusters with clinical labels and ",
          "Kaplan–Meier survival curves. PCA / UMAP coordinates are ",
          "precomputed; clustering re-runs on every change.")),
    layout_sidebar(
      sidebar = sidebar(
        width = 320, open = "always",
        h6("Selection"),
        selectInput(ns("cohort"),   "Cohort",   choices = NULL,
                    selected = DEFAULT_COHORT),
        selectInput(ns("modality"), "Modality", choices = NULL,
                    selected = DEFAULT_MODALITY),
        radioButtons(ns("tumor_only"), "Sample scope",
                     choices = c("Tumor only" = "tumor",
                                 "Tumor + normal" = "all"),
                     selected = "tumor", inline = TRUE),
        hr(),
        h6("Clustering"),
        radioButtons(ns("method"), "Method",
                     choices = c(`hierarchical (Ward.D2)` = "hierarchical",
                                 `k-means` = "kmeans"),
                     selected = "hierarchical"),
        sliderInput(ns("k"), "Number of clusters",
                    min = 2, max = 8, value = 3, step = 1),
        checkboxInput(ns("scale"), "Scale features (z-score)", TRUE),
        hr(),
        h6("Embedding"),
        radioButtons(ns("embed"), "Method",
                     choices = c("PCA" = "PCA", "UMAP" = "UMAP"),
                     selected = "PCA", inline = TRUE),
        selectInput(ns("color_by"), "Color points by",
                    choices = c("cluster","sample_type","stage","sex","OS_event"),
                    selected = "cluster"),
        hr(),
        downloadButton(ns("dl_clusters"), "Download cluster assignments (CSV)",
                       class = "btn-sm btn-outline-secondary"),
        downloadButton(ns("dl_km"),       "Save Kaplan–Meier plot (PNG)",
                       class = "btn-sm btn-outline-secondary")
      ),

      layout_columns(
        col_widths = c(8, 4),

        # ----------------------------------------------------- center -----
        tagList(
          layout_columns(
            col_widths = c(8, 4),
            card(card_header("Sample embedding"),
                 card_body(plotOutput(ns("embed_plot"), height = "440px"))),
            card(card_header("Cluster sizes"),
                 card_body(plotOutput(ns("size_plot"), height = "200px")),
                 card_body(plotOutput(ns("clinical_bar"), height = "200px")))
          ),
          layout_columns(
            col_widths = c(7, 5),
            card(card_header("Kaplan–Meier by cluster"),
                 card_body(plotOutput(ns("km_plot"), height = "360px"))),
            card(card_header("Top features by cluster"),
                 card_body(DT_or_table_ui(ns("markers_table"))))
          ),
          card(card_header("Cluster × marker heatmap"),
               card_body(plotOutput(ns("marker_heatmap"), height = "520px")))
        ),

        # ----------------------------------------------------- right rail -
        tagList(
          card(class = "right-rail",
               card_header("Interpretation"),
               card_body(uiOutput(ns("interp")))),
          card(class = "right-rail",
               card_header("Current solution"),
               card_body(uiOutput(ns("solution_summary")))),
          card(class = "right-rail",
               card_header("Reading the panels"),
               card_body(
                 p("The embedding reflects ", em("global"),
                   " sample geometry; the heatmap shows the cluster-defining ",
                   "features at z-score scale."),
                 p("KM curves include the log-rank p annotation. Cohorts ",
                   "with < 10 events have very wide confidence intervals."),
                 p(class = "text-caption",
                   "Cluster labels are unsupervised; treat them as a ",
                   "hypothesis to validate with stage/molecular subtype.")))
        )
      )
    )
  )
}

mod_molecular_structure_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {

    observe({
      m <- app_data$cohort_metadata
      choices <- if (!is.null(m))
        as.character(m[status == "completed", cohort]) else character()
      updateSelectInput(session, "cohort", choices = choices,
                        selected = isolate(input$cohort) %||% DEFAULT_COHORT)
    })

    observe({
      co <- input$cohort; req(co)
      mods <- modalities_with_cache(co)
      if (!length(mods)) mods <- names(MODALITIES)
      updateSelectInput(session, "modality",
                        choices = setNames(mods, modality_label(mods)),
                        selected = isolate(input$modality) %||% mods[1])
    })

    modality_obj <- reactive({
      req(input$cohort, input$modality)
      load_modality_cache(input$cohort, input$modality)
    })

    matrix_filtered <- reactive({
      mo <- modality_obj(); if (is.null(mo)) return(NULL)
      m  <- mo$matrix; ann <- as.data.table(mo$sample_meta)
      keep <- if (identical(input$tumor_only, "tumor"))
        ann$sample_type == "tumor" else rep(TRUE, nrow(ann))
      keep[is.na(keep)] <- FALSE
      list(matrix = m[, keep, drop = FALSE], ann = ann[keep])
    })

    clusters <- reactive({
      mf <- matrix_filtered()
      if (is.null(mf) || ncol(mf$matrix) < 4L) return(NULL)
      withProgress(message = sprintf("Clustering (%s, k=%d)…",
                                     input$method, input$k),
                   value = NULL, {
        run_clustering(mf$matrix,
                       method = input$method,
                       k      = input$k,
                       scale_features = input$scale)
      })
    })

    embedding <- reactive({
      mf <- matrix_filtered(); if (is.null(mf)) return(NULL)
      if (identical(input$embed, "PCA")) {
        pca <- load_pca_cache(input$cohort, input$modality)
        if (!is.null(pca) && !is.null(pca$pca_scores)) {
          scores <- pca$pca_scores
          df <- data.table(sample_id = rownames(scores),
                           dim1 = scores[, 1L], dim2 = scores[, 2L])
          df <- df[sample_id %in% colnames(mf$matrix)]
          return(df)
        }
      } else if (identical(input$embed, "UMAP")) {
        umap <- load_umap_cache(input$cohort, input$modality)
        if (!is.null(umap) && !is.null(umap$umap_scores)) {
          scores <- umap$umap_scores
          df <- data.table(sample_id = rownames(scores),
                           dim1 = scores[, 1L], dim2 = scores[, 2L])
          df <- df[sample_id %in% colnames(mf$matrix)]
          return(df)
        }
      }
      withProgress(message = sprintf("Computing %s embedding…", input$embed),
                   value = NULL, {
        embed_samples(mf$matrix, method = input$embed)
      })
    })

    embedding_with_meta <- reactive({
      e <- embedding(); cl <- clusters()
      mf <- matrix_filtered()
      if (is.null(e) || is.null(mf)) return(NULL)
      ann <- mf$ann
      e <- merge(e, ann[, intersect(c("sample_id","sample_type","stage",
                                       "sex","OS_event","OS_time","patient_id"),
                                     names(ann)), with = FALSE],
                 by = "sample_id", all.x = TRUE, sort = FALSE)
      if (!is.null(cl))
        e[, cluster := factor(unname(cl[match(sample_id, names(cl))]))]
      if ("OS_event" %in% names(e)) e[, OS_event := as.character(OS_event)]
      e
    })

    output$embed_plot <- renderPlot({
      d <- embedding_with_meta()
      if (is.null(d) || !nrow(d)) return(invisible(NULL))
      cb <- input$color_by
      if (!cb %in% names(d)) cb <- "cluster"
      title <- sprintf("%s · %s · %s · k=%d",
                       input$cohort,
                       modality_label(input$modality),
                       input$embed, input$k)
      p <- plot_embedding(d, color_by = cb, title = title)
      print(p)
    })

    output$size_plot   <- renderPlot({ print(plot_cluster_sizes(clusters())) })

    output$clinical_bar <- renderPlot({
      d <- embedding_with_meta(); if (is.null(d)) return(invisible(NULL))
      var_col <- if ("stage" %in% names(d) && any(!is.na(d$stage))) "stage"
                 else if ("sample_type" %in% names(d)) "sample_type"
                 else NULL
      if (is.null(var_col)) return(invisible(NULL))
      print(plot_cluster_vs_clinical(d, "cluster", var_col))
    })

    km_plot_obj <- reactive({
      d <- embedding_with_meta()
      if (is.null(d) || !all(c("OS_time","OS_event","cluster") %in% names(d)))
        return(NULL)
      df <- as.data.table(d)
      df[, OS_event := suppressWarnings(as.integer(OS_event))]
      plot_km_by_cluster(df)
    })
    output$km_plot <- renderPlot({
      p <- km_plot_obj(); if (is.null(p)) return(invisible(NULL)); print(p)
    })
    output$dl_km <- download_plot_handler(
      stem = function() sprintf("km_%s_%s_k%d",
                                  input$cohort %||% "cohort",
                                  input$modality %||% "modality",
                                  as.integer(input$k %||% 0L)),
      plot_reactive = km_plot_obj
    )

    markers <- reactive({
      mf <- matrix_filtered(); cl <- clusters()
      if (is.null(mf) || is.null(cl)) return(NULL)
      cluster_marker_table(mf$matrix, cl, top_n = 12L)
    })

    if (isTRUE(OPT_PKGS$DT)) {
      output$markers_table <- DT::renderDT({
        m <- markers(); if (is.null(m) || !nrow(m)) return(NULL)
        DT::datatable(m, rownames = FALSE,
                      options = list(pageLength = 10, dom = "ftip",
                                     scrollX = TRUE),
                      class = "compact stripe") |>
          DT::formatRound("diff_mean", digits = 3)
      })
    } else {
      output$markers_table <- renderTable({
        m <- markers(); if (is.null(m)) return(NULL); head(m, 30L)
      })
    }

    output$marker_heatmap <- renderPlot({
      mf <- matrix_filtered(); cl <- clusters(); mk <- markers()
      if (is.null(mf) || is.null(cl) || is.null(mk) || !nrow(mk))
        return(invisible(NULL))
      sel <- unique(head(mk[order(cluster, -diff_mean)], 80L)$feature)
      ann <- as.data.table(mf$ann)
      ann[, cluster := factor(unname(cl[match(sample_id, names(cl))]))]
      ann_cols <- intersect(c("cluster","sample_type","stage","sex"),
                             names(ann))
      m_sub <- mf$matrix[intersect(sel, rownames(mf$matrix)),
                          intersect(names(cl), colnames(mf$matrix)),
                          drop = FALSE]
      if (nrow(m_sub) < 2 || ncol(m_sub) < 4) return(invisible(NULL))
      plot_feature_heatmap(m_sub, sample_meta = ann,
                           ann_cols = ann_cols,
                           title = sprintf("%s · %s · top %d cluster markers",
                                           input$cohort,
                                           modality_label(input$modality),
                                           nrow(m_sub)))
    })

    # ---- right rail ----------------------------------------------------
    output$interp <- renderUI({
      mf <- matrix_filtered(); cl <- clusters()
      if (is.null(mf))
        return(empty_state("No modality cache",
                           "Run the data-prep script for this cohort/modality.",
                           icon = "exclamation-circle"))
      if (is.null(cl))
        return(empty_state("Insufficient samples",
                           "Not enough samples to cluster — try Tumor + normal scope.",
                           icon = "info-circle"))
      pca <- load_pca_cache(input$cohort, input$modality)
      pcv <- if (!is.null(pca)) pca$var_explained else NULL
      msgs <- interp_clustering(cl, mf$ann, pca_var_explained = pcv)
      tags$ul(class = "interp-list",
              lapply(msgs, function(m) tags$li(m)))
    })

    output$solution_summary <- renderUI({
      d <- embedding_with_meta(); cl <- clusters()
      if (is.null(d) || is.null(cl)) return(empty_state("Choose inputs",
                                                         "Pick a cohort/modality and clustering parameters."))
      tags$dl(class = "summary-dl",
        tags$dt("Cohort"),    tags$dd(input$cohort),
        tags$dt("Modality"),  tags$dd(modality_label(input$modality)),
        tags$dt("Scope"),     tags$dd(input$tumor_only),
        tags$dt("Method"),    tags$dd(input$method),
        tags$dt("k"),         tags$dd(input$k),
        tags$dt("Embedding"), tags$dd(input$embed),
        tags$dt("Color by"),  tags$dd(input$color_by),
        tags$dt("n samples"), tags$dd(length(cl)),
        tags$dt("Smallest cluster"), tags$dd(min(table(cl)))
      )
    })

    output$dl_clusters <- downloadHandler(
      filename = function() sprintf("%s_%s_clusters.csv",
                                    input$cohort, input$modality),
      content  = function(file) {
        d <- embedding_with_meta(); if (is.null(d)) d <- data.table()
        data.table::fwrite(d, file)
      }
    )
  })
}
