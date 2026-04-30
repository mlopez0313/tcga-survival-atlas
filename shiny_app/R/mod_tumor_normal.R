# mod_tumor_normal.R ----------------------------------------------------------
# Page 3 — Tumor vs Normal molecular biology.
# 3-pane layout: left controls / center evidence / right interpretation.
# -----------------------------------------------------------------------------

mod_tumor_normal_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "page-header",
        div(class = "page-eyebrow", "Modality view · per cohort"),
        h2("Tumor vs normal biology", class = "page-title"),
        p(class = "page-lede",
          "Compare RNA, miRNA, methylation, CNV (and mutation frequencies) ",
          "between tumor and matched-normal samples within a cohort. ",
          "All p-values are Welch's t-test on the top-variable feature ",
          "subset; FDR is BH-corrected.")),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, open = "always",
        h6("Selection"),
        selectInput(ns("cohort"), "Cohort",
                    choices = NULL, selected = DEFAULT_COHORT),
        selectInput(ns("modality"), "Modality",
                    choices = MODALITIES, selected = DEFAULT_MODALITY),
        selectizeInput(ns("feature"), "Feature",
                       choices = NULL, selected = NULL,
                       options = list(maxOptions = 200)),
        selectInput(ns("transform"), "Transform",
                    choices = c("none","z-score","log2(x+1)"),
                    selected = "none"),
        hr(),
        h6("Ranking"),
        radioButtons(ns("ranking"), NULL,
                     choices = c("absolute mean diff" = "abs_diff",
                                 "tumor − normal"     = "diff",
                                 "normal − tumor"     = "neg_diff",
                                 "smallest p-value"   = "pvalue"),
                     selected = "abs_diff"),
        sliderInput(ns("top_n"), "Top features in heatmap",
                    min = 10, max = 100, value = 30, step = 5),
        hr(),
        downloadButton(ns("dl_summary"),
                       "Download tumor/normal table (CSV)",
                       class = "btn-sm btn-outline-secondary"),
        downloadButton(ns("dl_volcano"),
                       "Save volcano plot (PNG)",
                       class = "btn-sm btn-outline-secondary")
      ),

      # ----- main content ----------------------------------------------------
      layout_columns(
        col_widths = c(8, 4),

        # center column ------------------------------------------------------
        tagList(
          layout_columns(
            col_widths = c(7, 5),
            card(card_header(textOutput(ns("feat_title"), inline = TRUE)),
                 card_body(plotOutput(ns("feat_box"), height = "360px"))),
            card(card_header("Summary stats"),
                 card_body(uiOutput(ns("feat_stats"))))
          ),
          card(card_header("Top differential features"),
               card_body(DT_or_table_ui(ns("tn_table")))),
          layout_columns(
            col_widths = c(7, 5),
            card(card_header("Volcano: tumor vs normal"),
                 card_body(plotOutput(ns("volcano"), height = "440px"))),
            card(card_header("PCA: tumor vs normal"),
                 card_body(plotOutput(ns("tn_pca"), height = "440px")))
          ),
          card(card_header("Top features × samples heatmap"),
               card_body(plotOutput(ns("feat_heatmap"), height = "520px"))),
          conditionalPanel(
            sprintf("input['%s'] == 'mutation'", ns("modality")),
            card(card_header("Mutation frequency (tumor only)"),
                 card_body(plotOutput(ns("mut_freq"), height = "560px"))))
        ),

        # right rail ---------------------------------------------------------
        tagList(
          card(class = "right-rail",
               card_header("Interpretation"),
               card_body(uiOutput(ns("interp")))),
          card(class = "right-rail",
               card_header("Cohort × modality at a glance"),
               card_body(uiOutput(ns("at_a_glance")))),
          card(class = "right-rail",
               card_header("Reading the volcano"),
               card_body(
                 p("X-axis: mean-difference between tumor and normal on the ",
                   "(transformed) modality scale. Y-axis: −log10 of Welch's ",
                   "t-test p."),
                 p("Dashed lines mark the auto-selected fold-change ",
                   "threshold (80th percentile of |Δ| for the displayed ",
                   "feature set) and p < 0.01."),
                 p(class = "text-caption",
                   "Methylation β-values live in [0,1] so |Δ| is small in ",
                   "absolute terms but biologically meaningful at ≥0.1.")))
        )
      )
    )
  )
}

mod_tumor_normal_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {

    # ---- choices --------------------------------------------------------
    observe({
      m <- app_data$cohort_metadata
      choices <- if (!is.null(m)) {
        ord <- order(m$status != "completed", m$cohort)
        setNames(m$cohort[ord],
                 ifelse(m$status[ord] == "completed",
                        m$cohort[ord],
                        paste0(m$cohort[ord], "  (", m$status[ord], ")")))
      } else character()
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

    has_normal <- reactive({
      mo <- modality_obj()
      if (is.null(mo)) return(FALSE)
      sum(mo$sample_meta$sample_type == "normal", na.rm = TRUE) >= 2L
    })

    # ---- top differential table ----------------------------------------
    tn_table <- reactive({
      req(input$cohort, input$modality)
      tn <- load_tn_summary(input$cohort, input$modality)
      if (is.null(tn)) tn <- compute_tn_on_fly(modality_obj())
      if (is.null(tn)) return(NULL)
      tn <- as.data.table(tn)
      tn <- switch(input$ranking,
                   abs_diff = tn[order(-abs(diff_mean))],
                   diff     = tn[order(-diff_mean)],
                   neg_diff = tn[order(diff_mean)],
                   pvalue   = if ("pvalue" %in% names(tn))
                                tn[order(pvalue, na.last = TRUE)]
                              else tn[order(-abs(diff_mean))])
      tn
    })

    # ---- feature picker -------------------------------------------------
    observe({
      tn <- tn_table()
      mo <- modality_obj()
      choices <- if (!is.null(tn) && nrow(tn)) tn$feature
                 else if (!is.null(mo)) rownames(mo$matrix)
                 else character()
      sel <- if (!is.null(isolate(input$feature)) &&
                 isolate(input$feature) %in% choices)
        isolate(input$feature) else (if (length(choices)) choices[1] else NULL)
      updateSelectizeInput(session, "feature",
                           choices = choices, selected = sel,
                           server = TRUE)
    })

    # ---- single-feature view -------------------------------------------
    feat_long <- reactive({
      mo <- modality_obj()
      if (is.null(mo) || !length(input$feature)) return(NULL)
      df <- feature_long_frame(mo, input$feature)
      if (is.null(df)) return(NULL)
      switch(input$transform,
             "z-score" = { v <- df$value
                           df[, value := (v - mean(v, na.rm = TRUE)) /
                                          sd(v, na.rm = TRUE)] },
             "log2(x+1)" = df[, value := log2(pmax(value, 0) + 1)],
             df)
      df
    })

    output$feat_title <- renderText({
      sprintf("%s — %s · %s",
              input$feature %||% "(no feature)",
              input$cohort %||% "—",
              modality_label(input$modality %||% ""))
    })

    output$feat_box <- renderPlot({
      df <- feat_long()
      if (is.null(df) || !nrow(df)) {
        return(invisible(NULL))
      }
      print(plot_feature_box(df,
                             feature_id = input$feature,
                             modality   = modality_label(input$modality)))
    })

    output$feat_stats <- renderUI({
      tn <- tn_table(); ft <- input$feature
      if (is.null(tn) || !length(ft))
        return(empty_state("Pick a feature",
                           "Search above or click any row in the top-features table."))
      row <- tn[feature == ft]
      if (!nrow(row))
        return(empty_state("No summary",
                           "This feature has no tumor/normal summary stats — likely all-NA in normals."))
      mo <- modality_obj()
      m  <- if (!is.null(mo)) mo$matrix[ft, , drop = TRUE] else NULL
      missing_pct <- if (!is.null(m)) 100 * mean(is.na(m)) else NA_real_
      pval_str <- if ("pvalue" %in% names(row) && !is.na(row$pvalue))
        sprintf("%.3g", row$pvalue) else "—"
      fdr_str <- if ("fdr" %in% names(row) && !is.na(row$fdr))
        sprintf("%.3g", row$fdr) else "—"
      tagList(
        tags$dl(class = "summary-dl",
          tags$dt("n tumor"),     tags$dd(row$n_tumor),
          tags$dt("n normal"),    tags$dd(row$n_normal),
          tags$dt("mean tumor"),  tags$dd(sprintf("%.3f", row$mean_tumor)),
          tags$dt("mean normal"), tags$dd(sprintf("%.3f", row$mean_normal)),
          tags$dt("median tumor"),  tags$dd(sprintf("%.3f", row$median_tumor)),
          tags$dt("median normal"), tags$dd(sprintf("%.3f", row$median_normal)),
          tags$dt("Δ mean (T−N)"),  tags$dd(sprintf("%.3f", row$diff_mean)),
          tags$dt("Welch p"),       tags$dd(pval_str),
          tags$dt("FDR"),           tags$dd(fdr_str),
          tags$dt("missingness"),   tags$dd(if (is.na(missing_pct)) "—"
                                            else sprintf("%.1f%%", missing_pct))
        )
      )
    })

    # ---- top differential table -----------------------------------------
    if (isTRUE(OPT_PKGS$DT)) {
      output$tn_table <- DT::renderDT({
        tn <- tn_table()
        if (is.null(tn) || !nrow(tn)) return(NULL)
        out <- head(tn, input$top_n %||% 30L)
        cols <- intersect(c("feature","mean_tumor","mean_normal","diff_mean",
                            "median_tumor","median_normal",
                            "pvalue","fdr",
                            "n_tumor","n_normal"), names(out))
        out <- out[, ..cols]
        dt <- DT::datatable(out, rownames = FALSE,
                      selection = list(mode = "single"),
                      options = list(pageLength = 10, dom = "ftip",
                                     scrollX = TRUE),
                      class = "compact stripe")
        round_cols <- intersect(c("mean_tumor","mean_normal","diff_mean",
                                   "median_tumor","median_normal"), names(out))
        if (length(round_cols)) dt <- DT::formatRound(dt, round_cols, digits = 3)
        sig_cols <- intersect(c("pvalue","fdr"), names(out))
        if (length(sig_cols)) dt <- DT::formatSignif(dt, sig_cols, digits = 3)
        dt
      })
      observeEvent(input$tn_table_rows_selected, {
        tn <- tn_table(); if (is.null(tn)) return()
        sel <- input$tn_table_rows_selected
        if (length(sel)) updateSelectizeInput(session, "feature",
                                              selected = tn$feature[sel[1]],
                                              server = TRUE)
      })
    } else {
      output$tn_table <- renderTable({
        tn <- tn_table(); if (is.null(tn)) return(NULL)
        head(tn, input$top_n %||% 30L)
      })
    }

    # ---- volcano --------------------------------------------------------
    volcano_plot <- reactive({
      tn <- tn_table()
      if (is.null(tn) || !nrow(tn) || !"pvalue" %in% names(tn)) return(NULL)
      fc_thresh <- max(0.5, quantile(abs(tn$diff_mean), 0.80, na.rm = TRUE))
      plot_volcano(tn, fc_thresh = fc_thresh, p_thresh = 0.01)
    })
    output$volcano <- renderPlot({
      p <- volcano_plot(); if (is.null(p)) return(invisible(NULL)); print(p)
    })
    output$dl_volcano <- download_plot_handler(
      stem = function() sprintf("volcano_%s_%s",
                                  input$cohort %||% "cohort",
                                  input$modality %||% "modality"),
      plot_reactive = volcano_plot
    )

    # ---- PCA ------------------------------------------------------------
    output$tn_pca <- renderPlot({
      mo <- modality_obj(); if (is.null(mo)) return(invisible(NULL))
      pca <- load_pca_cache(input$cohort, input$modality)
      if (!is.null(pca) && !is.null(pca$pca_scores)) {
        scores <- pca$pca_scores
        df <- data.table(sample_id = rownames(scores),
                         dim1 = scores[, 1L], dim2 = scores[, 2L])
        df[, sample_type := mo$sample_meta$sample_type[
          match(df$sample_id, mo$sample_meta$sample_id)]]
        v1 <- if (!is.null(pca$var_explained)) pca$var_explained[1] else NA
        v2 <- if (!is.null(pca$var_explained)) pca$var_explained[2] else NA
        title <- sprintf("PC1 (%.1f%%)  ·  PC2 (%.1f%%)", 100*v1, 100*v2)
      } else {
        df <- embed_samples(mo$matrix, method = "PCA")
        if (is.null(df)) return(invisible(NULL))
        df[, sample_type := mo$sample_meta$sample_type[
          match(df$sample_id, mo$sample_meta$sample_id)]]
        title <- "PCA (computed live)"
      }
      print(plot_embedding(df, color_by = "sample_type", title = title))
    })

    # ---- top features heatmap -------------------------------------------
    output$feat_heatmap <- renderPlot({
      mo <- modality_obj(); tn <- tn_table()
      if (is.null(mo) || is.null(tn) || !nrow(tn)) return(invisible(NULL))
      n_feat <- min(input$top_n %||% 30L, nrow(tn))
      feats <- head(tn, n_feat)$feature
      m_sub <- mo$matrix[intersect(feats, rownames(mo$matrix)), , drop = FALSE]
      ann <- as.data.table(mo$sample_meta)
      ann_cols <- intersect(c("sample_type","stage","sex"), names(ann))
      plot_feature_heatmap(m_sub, sample_meta = ann,
                           ann_cols = ann_cols,
                           title = sprintf("%s · %s · top %d differential",
                                           input$cohort,
                                           modality_label(input$modality),
                                           nrow(m_sub)))
    })

    # ---- mutation frequency --------------------------------------------
    output$mut_freq <- renderPlot({
      if (!identical(input$modality, "mutation")) return(invisible(NULL))
      mo <- modality_obj(); if (is.null(mo)) return(invisible(NULL))
      print(plot_mutation_freq(mo))
    })

    # ---- right rail: dynamic interpretation ----------------------------
    output$interp <- renderUI({
      mo <- modality_obj(); tn <- tn_table()
      if (is.null(mo))
        return(empty_state("No modality cache",
                           "Choose a different cohort or run the prep script for this modality.",
                           icon = "exclamation-circle"))
      msgs <- interp_tumor_normal(mo, tn, modality = input$modality)
      tags$ul(class = "interp-list",
              lapply(msgs, function(m) tags$li(m$children)))
    })

    output$at_a_glance <- renderUI({
      mo <- modality_obj()
      if (is.null(mo)) return(empty_state("No modality cache",
                                          "—", icon = "info-circle"))
      ann <- as.data.table(mo$sample_meta)
      n_t <- sum(ann$sample_type == "tumor",  na.rm = TRUE)
      n_n <- sum(ann$sample_type == "normal", na.rm = TRUE)
      tn <- tn_table()
      n_sig <- if (!is.null(tn) && "fdr" %in% names(tn))
        sum(tn$fdr < 0.01 & abs(tn$diff_mean) >= 1, na.rm = TRUE) else NA_integer_
      tags$dl(class = "summary-dl",
        tags$dt("Cohort"),   tags$dd(input$cohort),
        tags$dt("Modality"), tags$dd(modality_label(input$modality)),
        tags$dt("Features"), tags$dd(nrow(mo$matrix)),
        tags$dt("Samples"),  tags$dd(ncol(mo$matrix)),
        tags$dt("Tumor"),    tags$dd(n_t),
        tags$dt("Normal"),   tags$dd(n_n),
        tags$dt("FDR<1% & |Δ|≥1"), tags$dd(if (is.na(n_sig)) "—" else n_sig)
      )
    })

    output$dl_summary <- downloadHandler(
      filename = function() sprintf("%s_%s_tumor_normal.csv",
                                    input$cohort, input$modality),
      content  = function(file) {
        tn <- tn_table(); if (is.null(tn)) tn <- data.table()
        data.table::fwrite(tn, file)
      }
    )
  })
}
