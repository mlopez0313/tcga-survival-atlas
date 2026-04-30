# utils_plots.R ---------------------------------------------------------------
# ggplot2 / plotly helpers tailored for the app's editorial look.
# -----------------------------------------------------------------------------

#' Editorial ggplot theme.
theme_editorial <- function(base_size = 13) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major   = ggplot2::element_line(color = "#e3e8ee"),
      strip.text         = ggplot2::element_text(face = "bold"),
      plot.title         = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle      = ggplot2::element_text(color = "#5a6776"),
      axis.title         = ggplot2::element_text(color = "#2c3e50"),
      legend.position    = "right",
      plot.background    = ggplot2::element_rect(fill = "white", color = NA)
    )
}

#' Save a ggplot to a temp file in PNG (preferred) or PDF. Wired up by
#' `download_plot_handler()` below — no UI logic here.
save_ggplot <- function(plot, file, width = 9, height = 6, dpi = 200) {
  if (is.null(plot)) {
    grDevices::png(file, width = 600, height = 200, res = 100)
    on.exit(grDevices::dev.off())
    graphics::par(mar = c(2, 2, 2, 2)); graphics::plot.new()
    graphics::text(0.5, 0.5, "no plot available", cex = 1.2,
                   col = "#5a6776"); return(invisible())
  }
  ggplot2::ggsave(file, plot = plot, width = width, height = height,
                   dpi = dpi, bg = "white")
}

#' Build a downloadHandler for a reactive ggplot. `stem` may be a literal
#' string or a zero-arg function (closure) that returns one — the latter
#' lets the filename pick up live Shiny inputs at click time.
#' Usage:
#'    output$dl_plot <- download_plot_handler(
#'      stem = function() sprintf("volcano_%s_%s", input$cohort, input$modality),
#'      plot_reactive = volcano_plot)
download_plot_handler <- function(stem, plot_reactive,
                                   width = 9, height = 6) {
  shiny::downloadHandler(
    filename = function() {
      s <- if (is.function(stem)) stem() else stem
      if (!length(s) || !nzchar(s)) s <- "plot"
      sprintf("%s_%s.png", s, format(Sys.time(), "%Y%m%d-%H%M%S"))
    },
    contentType = "image/png",
    content = function(file) {
      p <- tryCatch(plot_reactive(), error = function(e) NULL)
      save_ggplot(p, file, width = width, height = height)
    })
}

#' Categorical palette used across the app (sourced from TCGA_PALETTE$qual).
app_palette <- function(n = 8) {
  base <- if (exists("TCGA_PALETTE")) TCGA_PALETTE$qual else
    c("#2c7fb8","#41ab5d","#fd8d3c","#756bb1","#fa9fb5",
      "#1d91c0","#d94801","#54278f","#bd0026","#525252")
  rep_len(base, n)
}

#' Consistent fills for tumor vs normal across the entire app.
tn_fill_scale <- function(...) {
  ggplot2::scale_fill_manual(
    values = c(tumor   = TCGA_PALETTE$tumor,
               normal  = TCGA_PALETTE$normal,
               control = TCGA_PALETTE$control),
    drop = FALSE, ...)
}
tn_color_scale <- function(...) {
  ggplot2::scale_color_manual(
    values = c(tumor   = TCGA_PALETTE$tumor,
               normal  = TCGA_PALETTE$normal,
               control = TCGA_PALETTE$control),
    drop = FALSE, ...)
}

#' Consistent colors for the four model families.
model_color_scale <- function(...) {
  vals <- TCGA_PALETTE$model_family
  names(vals) <- model_label(names(vals))
  ggplot2::scale_color_manual(values = vals, name = "Model", ...)
}
model_fill_scale <- function(...) {
  vals <- TCGA_PALETTE$model_family
  names(vals) <- model_label(names(vals))
  ggplot2::scale_fill_manual(values = vals, name = "Model", ...)
}

#' Wrap a ggplot in plotly if available, else return ggplot.
maybe_plotly <- function(p, tooltip = "all") {
  if (isTRUE(OPT_PKGS$plotly)) plotly::ggplotly(p, tooltip = tooltip) else p
}

#' Cohort × model heatmap of test C-index (numeric matrix expected).
plot_benchmark_heatmap <- function(bench_long) {
  if (is.null(bench_long) || !nrow(bench_long)) return(NULL)
  df <- as.data.table(bench_long)
  df[, model_label := model_label(model)]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = model_label, y = cohort,
                                         fill = c_index_test)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", c_index_test)),
                       size = 3, color = "#1f2d3d") +
    ggplot2::scale_fill_gradient2(low = "#deebf7", mid = "#9ecae1",
                                  high = "#08519c", midpoint = 0.65,
                                  limits = c(0.4, 1.0),
                                  name = "test C-index", oob = scales::squish) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Cross-cohort test C-index",
                  subtitle = "Higher = better discrimination") +
    theme_editorial() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  p
}

#' Win-count bar plot (uses the canonical per-model-family palette).
plot_win_counts <- function(win_counts) {
  if (is.null(win_counts) || !nrow(win_counts)) return(NULL)
  df <- as.data.table(win_counts)
  df[, model_label := model_label(model)]
  df[, model_label := factor(model_label, levels = model_label[order(-wins)])]
  ggplot2::ggplot(df, ggplot2::aes(x = wins, y = model_label,
                                    fill = model_label)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = wins), hjust = -0.3, size = 3.5,
                       color = TCGA_PALETTE$ink_soft) +
    model_fill_scale(guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(x = "# cohorts where model is best",
                  y = NULL, title = "Win counts across cohorts") +
    theme_editorial()
}

#' Tumor-vs-normal box+jitter plot for a single feature.
plot_feature_box <- function(feat_df, feature_id = "(feature)",
                              modality = "(modality)") {
  if (is.null(feat_df) || !nrow(feat_df)) return(NULL)
  df <- feat_df[!is.na(sample_type) & sample_type %in% c("tumor","normal")]
  if (!nrow(df)) return(NULL)
  df[, sample_type := factor(sample_type, levels = c("tumor","normal"))]
  ggplot2::ggplot(df, ggplot2::aes(x = sample_type, y = value,
                                    fill = sample_type, color = sample_type)) +
    ggplot2::geom_violin(alpha = 0.25, linewidth = 0.3, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85,
                          linewidth = 0.4, fill = "white") +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.4, size = 0.8) +
    tn_fill_scale(guide = "none") + tn_color_scale(guide = "none") +
    ggplot2::labs(x = NULL,
                  y = sprintf("%s (%s)", feature_id, modality),
                  title = sprintf("%s — tumor vs normal", feature_id),
                  subtitle = sprintf("n_tumor = %d, n_normal = %d",
                                     sum(df$sample_type == "tumor"),
                                     sum(df$sample_type == "normal"))) +
    theme_editorial()
}

#' PCA / UMAP scatter colored by an annotation column.
plot_embedding <- function(embed_df, color_by = NULL, title = NULL) {
  if (is.null(embed_df) || !nrow(embed_df)) return(NULL)
  aes_kw <- if (!is.null(color_by) && color_by %in% names(embed_df))
    ggplot2::aes(x = dim1, y = dim2, color = .data[[color_by]])
  else ggplot2::aes(x = dim1, y = dim2)
  p <- ggplot2::ggplot(embed_df, aes_kw) +
    ggplot2::geom_point(size = 1.7, alpha = 0.85) +
    ggplot2::labs(x = "Dim 1", y = "Dim 2", title = title,
                  color = if (!is.null(color_by)) color_by else NULL) +
    theme_editorial()
  if (!is.null(color_by) && color_by %in% names(embed_df) &&
      (is.character(embed_df[[color_by]]) || is.factor(embed_df[[color_by]]))) {
    if (color_by == "sample_type") {
      p <- p + tn_color_scale()
    } else {
      n_lev <- length(unique(embed_df[[color_by]]))
      p <- p + ggplot2::scale_color_manual(values = app_palette(n_lev))
    }
  }
  p
}

#' Cluster size bar plot.
plot_cluster_sizes <- function(clusters) {
  if (is.null(clusters) || !length(clusters)) return(NULL)
  df <- as.data.table(table(cluster = clusters))
  ggplot2::ggplot(df, ggplot2::aes(x = factor(cluster), y = N,
                                    fill = factor(cluster))) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.4, size = 3.5) +
    ggplot2::scale_fill_manual(values = app_palette(length(unique(clusters))),
                               guide = "none") +
    ggplot2::labs(x = "Cluster", y = "n samples",
                  title = "Cluster sizes") +
    theme_editorial()
}

#' Stacked bar of cluster vs categorical clinical variable.
plot_cluster_vs_clinical <- function(df, cluster_col = "cluster", var_col) {
  if (is.null(df) || !nrow(df) || !var_col %in% names(df)) return(NULL)
  d <- as.data.table(df)
  d <- d[!is.na(d[[var_col]])]
  if (!nrow(d)) return(NULL)
  ggplot2::ggplot(d, ggplot2::aes(x = .data[[cluster_col]],
                                    fill = .data[[var_col]])) +
    ggplot2::geom_bar(position = "fill", width = 0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_fill_manual(values = app_palette(
      length(unique(as.character(d[[var_col]]))))) +
    ggplot2::labs(x = "Cluster", y = "% of samples", fill = var_col,
                  title = sprintf("Cluster vs %s", var_col)) +
    theme_editorial()
}

#' Kaplan-Meier curves by cluster (uses survival pkg if available). Adds a
#' log-rank p annotation in the bottom-left so the panel is self-contained.
plot_km_by_cluster <- function(df) {
  if (is.null(df) || !nrow(df) || !OPT_PKGS$survival) return(NULL)
  df <- as.data.table(df)[!is.na(OS_time) & !is.na(OS_event)]
  df[, OS_event := suppressWarnings(as.integer(OS_event))]
  if (nrow(df) < 4 || sum(df$OS_event, na.rm = TRUE) < 2) return(NULL)
  fit <- tryCatch(survival::survfit(survival::Surv(OS_time, OS_event) ~ cluster,
                                     data = df),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  st <- summary(fit, censored = TRUE)
  curve <- data.table(time = st$time, surv = st$surv,
                      n.risk = st$n.risk,
                      strata = sub("^cluster=", "", as.character(st$strata)))
  # Log-rank
  lr <- tryCatch(
    survival::survdiff(survival::Surv(OS_time, OS_event) ~ cluster, data = df),
    error = function(e) NULL)
  p_lab <- if (!is.null(lr)) {
    pv <- if (!is.null(lr$pvalue)) lr$pvalue
          else 1 - pchisq(lr$chisq, df = max(length(lr$n) - 1L, 1L))
    sprintf("logrank p = %s",
            if (pv < 1e-3) format(pv, digits = 2, scientific = TRUE)
            else sprintf("%.3f", pv))
  } else NULL

  n_lev <- length(unique(curve$strata))
  pal   <- TCGA_PALETTE$qual[seq_len(n_lev)]

  p <- ggplot2::ggplot(curve, ggplot2::aes(x = time, y = surv,
                                            color = strata)) +
    ggplot2::geom_step(linewidth = 1.0) +
    ggplot2::scale_color_manual(values = pal, name = "Cluster") +
    ggplot2::labs(x = "Days", y = "Survival probability",
                  title = "Kaplan–Meier by cluster") +
    ggplot2::scale_y_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    theme_editorial()
  if (!is.null(p_lab)) {
    p <- p + ggplot2::annotate("label", x = -Inf, y = 0.05,
                                hjust = -0.05, vjust = 0,
                                label = p_lab, size = 3.4,
                                color = TCGA_PALETTE$ink_soft,
                                fill  = "white",
                                label.size = 0.2,
                                label.r = grid::unit(0.5, "lines"))
  }
  p
}

#' Mutation frequency bar plot (top genes).
plot_mutation_freq <- function(mut_obj, top_n = 25L) {
  if (is.null(mut_obj)) return(NULL)
  m <- mut_obj$matrix
  freq <- rowSums(m > 0, na.rm = TRUE)
  pct <- 100 * freq / ncol(m)
  df <- data.table(gene = names(freq), pct = pct, n = freq)[order(-pct)]
  df <- head(df, top_n)
  df[, gene := factor(gene, levels = rev(df$gene))]
  ggplot2::ggplot(df, ggplot2::aes(x = pct, y = gene)) +
    ggplot2::geom_col(fill = "#5d8aa8", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", pct)),
                       hjust = -0.2, size = 3) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(x = "% of tumor samples mutated", y = NULL,
                  title = sprintf("Top mutated genes (n = %d samples)",
                                  ncol(m))) +
    theme_editorial()
}

#' Volcano plot: log2-FC-like vs −log10(p) with significance shading.
plot_volcano <- function(tn_summary, fc_thresh = 1, p_thresh = 0.01,
                          top_label = 12L) {
  if (is.null(tn_summary) || !nrow(tn_summary)) return(NULL)
  d <- as.data.table(tn_summary)
  if (!all(c("diff_mean","pvalue") %in% names(d))) return(NULL)
  d <- d[is.finite(diff_mean) & is.finite(pvalue)]
  d[, neglog10_p := -log10(pmax(pvalue, .Machine$double.xmin))]
  d[, sig := fifelse(pvalue < p_thresh & abs(diff_mean) >= fc_thresh,
                     fifelse(diff_mean > 0, "up in tumor", "up in normal"),
                     "ns")]
  d[, sig := factor(sig, levels = c("up in tumor","up in normal","ns"))]
  to_label <- d[sig != "ns"][order(-neglog10_p)][seq_len(min(top_label, .N))]
  p <- ggplot2::ggplot(d, ggplot2::aes(x = diff_mean, y = neglog10_p,
                                        color = sig)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.4) +
    ggplot2::geom_hline(yintercept = -log10(p_thresh), linetype = "dashed",
                        color = "#888") +
    ggplot2::geom_vline(xintercept = c(-fc_thresh, fc_thresh),
                        linetype = "dashed", color = "#888") +
    ggplot2::scale_color_manual(values = c(`up in tumor`  = TCGA_PALETTE$tumor,
                                            `up in normal` = TCGA_PALETTE$normal,
                                            ns             = "#bdc3c7"),
                                 name = NULL) +
    ggplot2::labs(x = "Mean(tumor) − Mean(normal)",
                  y = expression(-log[10]~p),
                  title = "Volcano: tumor vs normal",
                  subtitle = sprintf("dashed: |diff| ≥ %.2g, p < %.2g",
                                     fc_thresh, p_thresh)) +
    theme_editorial()
  if (nrow(to_label) && OPT_PKGS$ggrepel) {
    p <- p + ggrepel::geom_text_repel(data = to_label,
                                       ggplot2::aes(label = feature),
                                       size = 3, max.overlaps = 20,
                                       segment.color = "#7f8c8d",
                                       show.legend = FALSE)
  } else if (nrow(to_label)) {
    p <- p + ggplot2::geom_text(data = to_label,
                                 ggplot2::aes(label = feature),
                                 size = 3, vjust = -0.6,
                                 show.legend = FALSE)
  }
  p
}

#' Per-feature heatmap (features × samples) annotated by sample_type/cluster.
#' Uses ComplexHeatmap if available, falls back to pheatmap, then to a
#' base-R image. `m` should already be a small subset (<= ~50 features).
plot_feature_heatmap <- function(m, sample_meta = NULL, ann_cols = NULL,
                                  scale_rows = TRUE, fontsize_row = 8,
                                  title = NULL) {
  if (is.null(m) || nrow(m) < 2 || ncol(m) < 4) return(invisible(NULL))
  X <- m
  if (scale_rows) {
    rsd <- apply(X, 1L, sd, na.rm = TRUE)
    keep <- which(is.finite(rsd) & rsd > 0)
    X <- X[keep, , drop = FALSE]
    X <- t(scale(t(X)))
    X[!is.finite(X)] <- 0
  }
  if (OPT_PKGS$ComplexHeatmap) {
    col_fun <- if (OPT_PKGS$circlize)
      circlize::colorRamp2(c(-2, 0, 2), c("#2c7fb8","white","#c0392b"))
    else NULL
    top_ann <- NULL
    if (!is.null(sample_meta) && length(ann_cols)) {
      sm <- as.data.table(sample_meta)[
        match(colnames(X), sample_meta$sample_id), , drop = FALSE]
      ann_keep <- intersect(ann_cols, names(sm))
      ann_df <- if (length(ann_keep))
        as.data.frame(sm[, ..ann_keep]) else data.frame()
      if (ncol(ann_df) >= 1L) {
        ann_pal <- lapply(names(ann_df), function(nm) {
          v <- as.character(ann_df[[nm]])
          lvls <- unique(v); lvls <- lvls[!is.na(lvls)]
          setNames(app_palette(length(lvls)), lvls)
        })
        names(ann_pal) <- names(ann_df)
        top_ann <- ComplexHeatmap::HeatmapAnnotation(df = ann_df,
                                                     col = ann_pal,
                                                     simple_anno_size = grid::unit(3.5, "mm"))
      }
    }
    ht <- ComplexHeatmap::Heatmap(
      X,
      name              = "z-score",
      col               = col_fun,
      top_annotation    = top_ann,
      show_row_names    = nrow(X) <= 80,
      show_column_names = FALSE,
      cluster_columns   = TRUE,
      cluster_rows      = TRUE,
      row_names_gp      = grid::gpar(fontsize = fontsize_row),
      column_title      = title,
      use_raster        = TRUE,
      raster_quality    = 2L
    )
    ComplexHeatmap::draw(ht, merge_legend = TRUE,
                         heatmap_legend_side = "right",
                         annotation_legend_side = "right")
    return(invisible(NULL))
  }
  if (OPT_PKGS$pheatmap) {
    ann_df <- NULL
    if (!is.null(sample_meta) && length(ann_cols)) {
      sm <- as.data.table(sample_meta)[
        match(colnames(X), sample_meta$sample_id), , drop = FALSE]
      ann_keep <- intersect(ann_cols, names(sm))
      if (length(ann_keep)) {
        ann_df <- as.data.frame(sm[, ..ann_keep])
        rownames(ann_df) <- colnames(X)
      }
    }
    pheatmap::pheatmap(
      X,
      annotation_col = ann_df,
      show_rownames  = nrow(X) <= 80,
      show_colnames  = FALSE,
      fontsize_row   = fontsize_row,
      color = colorRampPalette(c("#2c7fb8","white","#c0392b"))(99),
      main  = title %||% ""
    )
    return(invisible(NULL))
  }
  graphics::image(t(X), col = grDevices::hcl.colors(64, "RdBu", rev = TRUE),
                  axes = FALSE, main = title)
}

#' Modality availability heatmap across cohorts.
plot_modality_heatmap <- function(meta) {
  if (is.null(meta)) return(NULL)
  cols <- c("modality_rna","modality_mirna","modality_methylation",
            "modality_cnv","modality_mutation","modality_clinical")
  if (!all(cols %in% names(meta))) return(NULL)
  long <- melt(as.data.table(meta)[, c("cohort","status", cols), with = FALSE],
               id.vars = c("cohort","status"),
               measure.vars = cols, variable.name = "modality",
               value.name = "present")
  long[, modality := sub("^modality_", "", modality)]
  long[, modality := factor(modality, levels = c("clinical","rna","mirna",
                                                 "methylation","cnv","mutation"))]
  long[, present := as.logical(present)]
  ggplot2::ggplot(long, ggplot2::aes(x = modality, y = cohort,
                                      fill = present)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#2c7fb8",
                                          `FALSE` = "#e3e8ee"),
                               labels = c(`TRUE` = "present",
                                          `FALSE` = "missing"),
                               name = NULL) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::labs(x = NULL, y = NULL, title = "Modality availability") +
    theme_editorial() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
}

#' Performance vs sample count scatter (model family colors are consistent
#' with the rest of the app via `model_color_scale()`).
plot_perf_vs_size <- function(bench_long) {
  if (is.null(bench_long) || !nrow(bench_long)) return(NULL)
  df <- as.data.table(bench_long)
  df[, model_label := model_label(model)]
  ggplot2::ggplot(df, ggplot2::aes(x = events_test + events_train,
                                    y = c_index_test, color = model_label,
                                    label = cohort)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed",
                        color = "#b0b8c1", linewidth = 0.4) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    model_color_scale() +
    ggplot2::labs(x = "Total events (train + test)",
                  y = "Test C-index",
                  title = "Performance vs event count",
                  subtitle = "0.5 dashed line = no discrimination") +
    theme_editorial()
}
