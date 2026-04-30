# utils_data_prep.R -----------------------------------------------------------
# Lightweight data shaping for the modules.
# -----------------------------------------------------------------------------

#' Null-coalescing operator. Used widely by the module servers to fall back to
#' sane defaults from selectInput choices that haven't initialized yet.
`%||%` <- function(a, b) {
  if (length(a) == 0L) return(b)
  if (length(a) == 1L) {
    if (is.na(a)) return(b)
    if (is.character(a) && !nzchar(a)) return(b)
  }
  a
}

#' Build a wide cohort × model matrix of test C-index from `benchmark_long`.
#' Returns a data.frame; rows = cohort, columns = model.
benchmark_wide <- function(benchmark_long) {
  if (is.null(benchmark_long) || !nrow(benchmark_long)) return(NULL)
  bl <- as.data.table(benchmark_long)
  if (!"model" %in% names(bl) || !"c_index_test" %in% names(bl)) return(NULL)
  dcast(bl, cohort ~ model, value.var = "c_index_test", fun.aggregate = mean)
}

#' Pretty model label (falls back to raw value).
model_label <- function(x) {
  unname(ifelse(x %in% names(MODEL_LABELS), MODEL_LABELS[x], x))
}

#' Pretty modality label.
modality_label <- function(x) {
  unname(ifelse(x %in% names(MODALITIES), MODALITIES[x], x))
}

#' Long-format frame for a single feature × samples (tumor/normal in long form).
#' Returns data.table with cols: sample_id, sample_type, value.
feature_long_frame <- function(modality_obj, feature_id) {
  if (is.null(modality_obj)) return(NULL)
  m <- modality_obj$matrix
  if (!feature_id %in% rownames(m)) return(NULL)
  vals <- m[feature_id, , drop = TRUE]
  data.table(
    sample_id   = colnames(m),
    sample_type = modality_obj$sample_meta$sample_type,
    patient_id  = modality_obj$sample_meta$patient_id,
    value       = as.numeric(vals)
  )
}

#' Compute tumor-vs-normal summary on the fly when no cache exists.
compute_tn_on_fly <- function(modality_obj, top_n = NULL) {
  if (is.null(modality_obj)) return(NULL)
  m  <- modality_obj$matrix
  st <- modality_obj$sample_meta$sample_type
  is_t <- st == "tumor"; is_n <- st == "normal"
  if (sum(is_t) < 2 || sum(is_n) < 2) return(NULL)

  rm <- function(x, idx) rowMeans(x[, idx, drop = FALSE], na.rm = TRUE)
  rmd <- function(x, idx) apply(x[, idx, drop = FALSE], 1L, median, na.rm = TRUE)
  rsd <- function(x, idx) apply(x[, idx, drop = FALSE], 1L, sd, na.rm = TRUE)

  out <- data.table(
    feature       = rownames(m),
    n_tumor       = sum(is_t),
    n_normal      = sum(is_n),
    mean_tumor    = rm(m, which(is_t)),
    mean_normal   = rm(m, which(is_n)),
    median_tumor  = rmd(m, which(is_t)),
    median_normal = rmd(m, which(is_n)),
    sd_tumor      = rsd(m, which(is_t)),
    sd_normal     = rsd(m, which(is_n))
  )
  out[, diff_mean := mean_tumor - mean_normal]
  out[, log2fc_like := diff_mean]
  out <- out[order(-abs(diff_mean))]
  if (!is.null(top_n)) head(out, top_n) else out
}

#' Run k-means or hierarchical clustering on a (features × samples) matrix.
#' Returns integer cluster labels of length ncol(m).
run_clustering <- function(m, method = c("hierarchical","kmeans"), k = 3L,
                            scale_features = TRUE, max_features = 2000L) {
  method <- match.arg(method)
  if (is.null(m) || ncol(m) < 3 || nrow(m) < 2) return(NULL)
  X <- t(m)                                      # samples × features
  X[!is.finite(X)] <- NA
  col_means <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) {
    nas <- is.na(X[, j]); if (any(nas)) X[nas, j] <- col_means[j]
  }
  if (scale_features) X <- scale(X)
  X[!is.finite(X)] <- 0
  if (ncol(X) > max_features) {
    v <- apply(X, 2L, var, na.rm = TRUE)
    X <- X[, order(v, decreasing = TRUE)[seq_len(max_features)], drop = FALSE]
  }
  if (method == "kmeans") {
    set.seed(42L)
    km <- tryCatch(kmeans(X, centers = k, nstart = 10L, iter.max = 50L),
                   error = function(e) NULL)
    if (is.null(km)) return(NULL)
    setNames(as.integer(km$cluster), rownames(X))
  } else {
    d <- tryCatch(dist(X), error = function(e) NULL)
    if (is.null(d)) return(NULL)
    hc <- hclust(d, method = "ward.D2")
    setNames(as.integer(cutree(hc, k = k)), rownames(X))
  }
}

#' Compute UMAP/PCA embedding from features × samples matrix.
#' Returns data.table(sample_id, dim1, dim2) or NULL.
embed_samples <- function(m, method = c("PCA","UMAP"),
                          n_neighbors = 15L, min_dist = 0.1) {
  method <- match.arg(method)
  if (is.null(m) || ncol(m) < 4) return(NULL)
  X <- t(m); X[!is.finite(X)] <- NA
  cm <- colMeans(X, na.rm = TRUE)
  for (j in seq_len(ncol(X))) {
    nas <- is.na(X[, j]); if (any(nas)) X[nas, j] <- cm[j]
  }
  v <- apply(X, 2L, sd, na.rm = TRUE)
  X <- X[, v > 0 & is.finite(v), drop = FALSE]
  if (ncol(X) < 2) return(NULL)
  if (method == "PCA") {
    pr <- tryCatch(prcomp(X, center = TRUE, scale. = TRUE,
                          rank. = min(2L, ncol(X), nrow(X) - 1L)),
                   error = function(e) NULL)
    if (is.null(pr)) return(NULL)
    data.table(sample_id = rownames(X),
               dim1 = pr$x[, 1L],
               dim2 = pr$x[, 2L])
  } else {
    if (!OPT_PKGS$uwot) return(NULL)
    set.seed(42L)
    u <- tryCatch(uwot::umap(scale(X),
                             n_neighbors = min(n_neighbors, nrow(X) - 1L),
                             min_dist    = min_dist,
                             n_components= 2L,
                             verbose     = FALSE),
                  error = function(e) NULL)
    if (is.null(u)) return(NULL)
    data.table(sample_id = rownames(X),
               dim1 = u[, 1L], dim2 = u[, 2L])
  }
}

#' Per-cluster KM data assembled from sample annotations + cluster labels.
#' Returns a list(time, event, cluster, df) or NULL when too few events.
km_data_for_clusters <- function(sample_ann, clusters) {
  if (is.null(sample_ann) || is.null(clusters)) return(NULL)
  if (!"OS_time" %in% names(sample_ann) || !"OS_event" %in% names(sample_ann))
    return(NULL)
  df <- as.data.table(sample_ann)
  df <- df[match(names(clusters), df$sample_id)]
  df[, cluster := factor(clusters)]
  df <- df[!is.na(OS_time) & !is.na(OS_event)]
  if (nrow(df) < 4 || sum(df$OS_event) < 2) return(NULL)
  list(df = df)
}

#' Calm empty-state card body. Keeps tone consistent across modules.
empty_state <- function(title, body = NULL, icon = "info-circle") {
  htmltools::div(
    class = "empty-state",
    htmltools::div(class = "empty-state-icon",
                   bsicons::bs_icon(icon)),
    htmltools::h6(title),
    if (!is.null(body)) htmltools::p(class = "empty-state-body", body)
    else NULL
  )
}

#' One-paragraph dynamic interpretation for the tumor-vs-normal view.
interp_tumor_normal <- function(modality_obj, tn_summary, modality = "") {
  if (is.null(modality_obj)) {
    return(list(p("No data loaded for this cohort/modality.")))
  }
  ann <- as.data.table(modality_obj$sample_meta)
  n_t <- sum(ann$sample_type == "tumor",  na.rm = TRUE)
  n_n <- sum(ann$sample_type == "normal", na.rm = TRUE)
  msgs <- list()
  msgs[[length(msgs) + 1L]] <-
    sprintf("This view shows %d tumor and %d normal samples.", n_t, n_n)
  if (n_n < 2L) {
    msgs[[length(msgs) + 1L]] <- paste0(
      "No matched-normal samples are available, so per-feature tumor-vs-",
      "normal comparisons are disabled. The PCA plot still shows the ",
      "structure of the tumor-only sample population.")
    return(lapply(msgs, p))
  }
  if (!is.null(tn_summary) && nrow(tn_summary) && "fdr" %in% names(tn_summary)) {
    n_sig <- sum(tn_summary$fdr < 0.01 & abs(tn_summary$diff_mean) >= 1, na.rm = TRUE)
    n_total <- nrow(tn_summary)
    pct <- 100 * n_sig / max(n_total, 1L)
    if (pct >= 30) {
      msgs[[length(msgs) + 1L]] <- sprintf(paste0(
        "Tumor–normal separation is visually strong: %d of %d top-variable ",
        "features (%.0f%%) cross the |Δ|≥1 / FDR<1%% threshold."),
        n_sig, n_total, pct)
    } else if (pct >= 5) {
      msgs[[length(msgs) + 1L]] <- sprintf(paste0(
        "Tumor–normal separation is moderate: %d of %d (%.0f%%) features ",
        "cross the |Δ|≥1 / FDR<1%% threshold."),
        n_sig, n_total, pct)
    } else {
      msgs[[length(msgs) + 1L]] <- sprintf(paste0(
        "Tumor–normal separation is subtle: only %d of %d (%.1f%%) features ",
        "cross the |Δ|≥1 / FDR<1%% threshold. Read individual hits with care."),
        n_sig, n_total, pct)
    }
  }
  if (modality == "mutation") {
    msgs[[length(msgs) + 1L]] <- paste0(
      "Mutation calls only exist in tumor samples, so the volcano / box ",
      "comparisons are intentionally disabled — see the frequency view below.")
  }
  if (modality == "methylation") {
    msgs[[length(msgs) + 1L]] <- paste0(
      "β-values are bounded in [0,1]; the y-axis on box plots and the ",
      "diff_mean column should be read as fractional methylation change.")
  }
  lapply(msgs, p)
}

#' Dynamic interpretation for a clustering solution.
interp_clustering <- function(clusters, sample_ann, embedding_df = NULL,
                               pca_var_explained = NULL) {
  if (is.null(clusters) || !length(clusters))
    return(list(p("Cluster solution not available.")))
  msgs <- list()
  k <- length(unique(clusters))
  smallest <- min(table(clusters))
  msgs[[length(msgs) + 1L]] <-
    sprintf("Solution: k=%d (smallest cluster: n=%d).", k, smallest)
  if (smallest < 5L) {
    msgs[[length(msgs) + 1L]] <- paste0(
      "At least one cluster has fewer than 5 samples — interpret as a ",
      "tentative substructure, not a robust subtype.")
  }
  if (!is.null(sample_ann)) {
    sa <- as.data.table(sample_ann)
    sa <- sa[match(names(clusters), sa$sample_id)]
    sa[, cluster := factor(unname(clusters))]
    # Stage alignment
    if ("stage" %in% names(sa) && sum(!is.na(sa$stage)) > 10L) {
      tab <- table(sa$cluster, sa$stage)
      if (all(dim(tab) >= 2L)) {
        chi <- suppressWarnings(chisq.test(tab))
        if (!is.na(chi$p.value)) {
          msg <- if (chi$p.value < 0.01)
            sprintf("Clusters track stage strongly (χ² p = %.2g).", chi$p.value)
          else if (chi$p.value < 0.05)
            sprintf("Clusters and stage are weakly associated (χ² p = %.2g).",
                    chi$p.value)
          else
            "Clusters do not align cleanly with stage."
          msgs[[length(msgs) + 1L]] <- msg
        }
      }
    }
    # Survival separation via logrank
    if (all(c("OS_time","OS_event") %in% names(sa)) && OPT_PKGS$survival) {
      df <- sa[!is.na(OS_time) & !is.na(OS_event)]
      df[, OS_event := suppressWarnings(as.integer(OS_event))]
      if (nrow(df) >= 20L && length(unique(df$cluster)) >= 2L &&
          sum(df$OS_event, na.rm = TRUE) >= 5L) {
        lr <- tryCatch(
          survival::survdiff(survival::Surv(OS_time, OS_event) ~ cluster,
                              data = df),
          error = function(e) NULL
        )
        if (!is.null(lr)) {
          p_lr <- if (!is.null(lr$pvalue)) lr$pvalue
                  else 1 - pchisq(lr$chisq, df = length(lr$n) - 1L)
          msg <- if (p_lr < 0.01)
            sprintf("Survival differs strongly between clusters (logrank p = %.2g).", p_lr)
          else if (p_lr < 0.05)
            sprintf("Survival differs between clusters (logrank p = %.2g).", p_lr)
          else
            "This cluster solution does not separate survival cleanly."
          msgs[[length(msgs) + 1L]] <- msg
          attr(msgs, "logrank_p") <- p_lr
        }
      }
    }
    # Tumor/normal mix
    if ("sample_type" %in% names(sa)) {
      n_t <- sum(sa$sample_type == "tumor",  na.rm = TRUE)
      n_n <- sum(sa$sample_type == "normal", na.rm = TRUE)
      if (n_n >= 2L && n_t >= 2L) {
        # Are normals concentrated in one cluster?
        norm_clust <- table(sa$cluster[sa$sample_type == "normal"])
        if (length(norm_clust) && max(norm_clust) / sum(norm_clust) >= 0.8) {
          dom <- names(norm_clust)[which.max(norm_clust)]
          msgs[[length(msgs) + 1L]] <- sprintf(
            "Normals concentrate in cluster %s (%d/%d) — clusters likely reflect tumor-vs-normal first.",
            dom, max(norm_clust), sum(norm_clust))
        }
      }
    }
  }
  if (!is.null(pca_var_explained)) {
    pc12 <- 100 * sum(pca_var_explained[1:2], na.rm = TRUE)
    if (pc12 >= 50)
      msgs[[length(msgs) + 1L]] <- sprintf(
        "PC1+PC2 capture %.0f%% of the variance — the 2D embedding is highly informative.",
        pc12)
  }
  msgs
}

#' Dynamic interpretation for a single cohort's survival benchmark.
interp_survival_cohort <- function(summary_df, cohort) {
  if (is.null(summary_df) || !nrow(summary_df))
    return(list(p("No model summary available.")))
  d <- as.data.table(summary_df)
  msgs <- list()
  best <- d[which.max(c_index_test)]
  msgs[[length(msgs) + 1L]] <- sprintf(
    "Best model for %s: %s (test C-index %.3f).",
    cohort, model_label(best$model[1]), best$c_index_test[1])
  spread <- diff(range(d$c_index_test, na.rm = TRUE))
  if (spread < 0.02) {
    msgs[[length(msgs) + 1L]] <- "Model performance is virtually identical (Δ C < 0.02); model choice is not meaningful here."
  } else if (spread < 0.05) {
    msgs[[length(msgs) + 1L]] <- "Models cluster within 0.05 C-index; differences are inside the noise floor for this cohort."
  }
  if ("events_test" %in% names(d) && !is.na(d$events_test[1]) &&
      d$events_test[1] < 30) {
    msgs[[length(msgs) + 1L]] <- sprintf(
      "This cohort has %d test events — interpret survival separations cautiously; CIs will be wide.",
      d$events_test[1])
  }
  msgs
}

#' Per-cluster top differential features (one-vs-rest mean difference).
cluster_marker_table <- function(m, clusters, top_n = 20L) {
  if (is.null(m) || is.null(clusters)) return(NULL)
  shared <- intersect(colnames(m), names(clusters))
  if (length(shared) < 4) return(NULL)
  m <- m[, shared, drop = FALSE]
  cl <- clusters[shared]
  out <- list()
  for (k in unique(cl)) {
    in_k  <- which(cl == k); out_k <- which(cl != k)
    if (length(in_k) < 2 || length(out_k) < 2) next
    diff <- rowMeans(m[, in_k, drop = FALSE], na.rm = TRUE) -
            rowMeans(m[, out_k, drop = FALSE], na.rm = TRUE)
    out[[as.character(k)]] <- data.table(
      feature      = names(diff),
      cluster      = k,
      diff_mean    = unname(diff)
    )[order(-diff_mean)][seq_len(min(top_n, length(diff)))]
  }
  if (!length(out)) return(NULL)
  rbindlist(out, use.names = TRUE, fill = TRUE)
}
