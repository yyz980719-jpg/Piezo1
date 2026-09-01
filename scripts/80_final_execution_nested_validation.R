# PIEZO1 OA final execution: truly nested leave-one-donor-out validation
# and high-recall-state sensitivity analysis. R-only analytical workflow.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

set.seed(20260831)
source("scripts/_bootstrap.R")
ROOT <- PIEZO1_ROOT
PRIOR <- PIEZO1_METADATA
OUT <- PIEZO1_METADATA
FIG <- PIEZO1_FIGURES
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

balanced_accuracy <- function(truth, pred, lev) {
  recalls <- vapply(lev, function(s) {
    ii <- truth == s
    if (!any(ii)) return(NA_real_)
    mean(pred[ii] == s)
  }, numeric(1))
  mean(recalls, na.rm = TRUE)
}

macro_f1 <- function(truth, pred, lev) {
  f1 <- vapply(lev, function(s) {
    tp <- sum(truth == s & pred == s)
    fp <- sum(truth != s & pred == s)
    fn <- sum(truth == s & pred != s)
    precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    recall <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  }, numeric(1))
  mean(f1)
}

message("Loading frozen reference object")
seu <- readRDS(file.path(ROOT, "data/scrna/GSE104782_seurat.rds"))
md <- as.data.table(seu@meta.data, keep.rownames = "cell")
md[, donor := as.character(patient)]
md[, truth := as.character(subtype)]
log_expr <- tryCatch(GetAssayData(seu, layer = "data"),
                     error = function(e) GetAssayData(seu, slot = "data"))
linear_expr <- log_expr
linear_expr@x <- expm1(linear_expr@x)
states <- sort(unique(md$truth))
donors <- sort(unique(md$donor))
genes <- rownames(log_expr)
n_cells <- ncol(log_expr)

# Query standardization is independent of the reference labels, so it is
# precomputed once per held-out donor. Marker selection and centroids are not.
query_z <- lapply(donors, function(dd) {
  ii <- which(md$donor == dd)
  q <- as.matrix(log_expr[, ii, drop = FALSE])
  z <- t(scale(t(q)))
  z[!is.finite(z)] <- 0
  z
})
names(query_z) <- donors

group_levels <- as.vector(outer(donors, states, paste, sep = "||"))

run_nested <- function(labels, keep_details = FALSE) {
  group_id <- match(paste(md$donor, labels, sep = "||"), group_levels)
  G <- sparseMatrix(i = seq_len(n_cells), j = group_id, x = 1,
                    dims = c(n_cells, length(group_levels)))
  sums <- as.matrix(linear_expr %*% G)
  counts <- tabulate(group_id, nbins = length(group_levels))
  colnames(sums) <- group_levels

  total_sums <- sapply(states, function(s) {
    rowSums(sums[, paste(donors, s, sep = "||"), drop = FALSE])
  })
  total_counts <- vapply(states, function(s) {
    sum(counts[match(paste(donors, s, sep = "||"), group_levels)])
  }, numeric(1))

  prediction_parts <- vector("list", length(donors))
  marker_parts <- vector("list", length(donors))
  for (k in seq_along(donors)) {
    dd <- donors[k]
    held_cols <- match(paste(dd, states, sep = "||"), group_levels)
    train_counts <- total_counts - counts[held_cols]
    train_means <- sweep(total_sums - sums[, held_cols, drop = FALSE],
                         2, train_counts, "/")
    colnames(train_means) <- states

    marker_rows <- vector("list", length(states))
    for (j in seq_along(states)) {
      other <- rowMeans(train_means[, -j, drop = FALSE])
      fc <- (train_means[, j] + 1e-6) / (other + 1e-6)
      fc[train_means[, j] < 0.05] <- 0
      ord <- order(fc, decreasing = TRUE)[seq_len(100)]
      marker_rows[[j]] <- data.table(fold_donor = dd, state = states[j],
                                     gene = genes[ord], fold_change = fc[ord])
    }
    fold_markers <- rbindlist(marker_rows)
    selected <- unique(fold_markers$gene)
    centroids <- train_means[match(selected, genes), , drop = FALSE]
    centroid_z <- t(scale(t(centroids)))
    centroid_z[!is.finite(centroid_z)] <- 0
    scores <- t(query_z[[dd]][match(selected, genes), , drop = FALSE]) %*% centroid_z
    pred <- states[max.col(scores, ties.method = "first")]
    ii <- which(md$donor == dd)
    prediction_parts[[k]] <- data.table(cell = md$cell[ii], donor = dd,
                                        truth = labels[ii], pred = pred)
    if (keep_details) marker_parts[[k]] <- fold_markers
  }
  predictions <- rbindlist(prediction_parts)
  result <- list(
    balanced_accuracy = balanced_accuracy(predictions$truth, predictions$pred, states),
    macro_f1 = macro_f1(predictions$truth, predictions$pred, states)
  )
  if (keep_details) {
    result$predictions <- predictions
    result$markers <- rbindlist(marker_parts)
  }
  result
}

message("Running observed fully nested LODO pipeline")
observed <- run_nested(md$truth, keep_details = TRUE)

# Pre-generate donor-preserving permutations on the master process. Each
# permutation repeats marker selection, centroid construction, prediction and
# scoring within every held-out donor fold.
n_perm <- 1000L
perm_labels <- lapply(seq_len(n_perm), function(i) {
  unlist(lapply(donors, function(dd) {
    ii <- which(md$donor == dd)
    sample(md$truth[ii], length(ii), replace = FALSE)
  }), use.names = FALSE)
})
# The donor blocks are already in donor-sorted order in this object; map them
# back to original cell order explicitly to guard against future metadata order.
donor_order <- unlist(lapply(donors, function(dd) which(md$donor == dd)))
perm_labels <- lapply(perm_labels, function(x) {
  y <- character(n_cells)
  y[donor_order] <- x
  y
})

workers <- max(1L, min(8L, parallel::detectCores(logical = FALSE) - 1L))
message(sprintf("Running %d complete donor-preserving permutations on %d workers", n_perm, workers))
cl <- parallel::makePSOCKcluster(workers)
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({library(Matrix); library(data.table)})
  NULL
})
parallel::clusterExport(
  cl,
  c("md", "log_expr", "linear_expr", "states", "donors", "genes", "n_cells",
    "query_z", "group_levels", "balanced_accuracy", "macro_f1", "run_nested"),
  envir = environment()
)
null_ba <- numeric(n_perm)
block_size <- 25L
for (lo in seq(1L, n_perm, by = block_size)) {
  hi <- min(n_perm, lo + block_size - 1L)
  ans <- parallel::parLapplyLB(cl, perm_labels[lo:hi], function(z) run_nested(z, FALSE))
  null_ba[lo:hi] <- vapply(ans, `[[`, numeric(1), "balanced_accuracy")
  message(sprintf("Completed permutations %d-%d", lo, hi))
}
parallel::stopCluster(cl)

pred <- observed$predictions
conf <- as.data.table(table(factor(pred$truth, levels = states),
                            factor(pred$pred, levels = states)))
setnames(conf, c("truth", "pred", "N"))
conf[, row_total := sum(N), by = truth]
conf[, row_fraction := N / row_total]
recall <- pred[, .(estimate = mean(pred == truth), n = .N), by = truth]
null95 <- unname(quantile(null_ba, 0.95, type = 8))
perm_p <- (1 + sum(null_ba >= observed$balanced_accuracy)) / (n_perm + 1)
qc <- rbind(
  data.table(metric = "balanced_accuracy", state = NA_character_,
             estimate = observed$balanced_accuracy, n = nrow(pred),
             null_95th = null95, permutation_p = perm_p,
             pass = observed$balanced_accuracy > null95),
  data.table(metric = "macro_F1", state = NA_character_,
             estimate = observed$macro_f1, n = nrow(pred),
             null_95th = NA_real_, permutation_p = NA_real_, pass = NA),
  recall[, .(metric = "recall", state = truth, estimate, n,
             null_95th = NA_real_, permutation_p = NA_real_, pass = estimate >= 0.60)]
)

fwrite(qc, file.path(OUT, "04_Nested_LODO_QC.csv"))
fwrite(conf, file.path(OUT, "04_LODO_confusion_matrix.csv"))
fwrite(pred, file.path(OUT, "04_Nested_LODO_CellPredictions.csv"))
fwrite(observed$markers, file.path(OUT, "04_Nested_LODO_Markers.csv"))
fwrite(data.table(iteration = seq_len(n_perm), balanced_accuracy = null_ba),
       file.path(OUT, "04_Nested_LODO_PermutationNull.csv"))

# High-recall-state-only sensitivity analysis. The transferred labels remain the
# frozen full-reference assignments; this step only excludes states that failed
# the nested recall threshold and renormalizes retained-state proportions.
assignment <- fread(file.path(PRIOR, "01_GSE152805_StateAssignmentConfidence.csv"))
low_recall <- recall[estimate < 0.60, truth]
retained <- setdiff(states, low_recall)
sens <- assignment[state %in% retained]

decomp_rows <- list()
direction_rows <- list()
for (dd in sort(unique(sens$donor))) {
  la <- sens[donor == dd & compartment == "lateral_nonload"]
  me <- sens[donor == dd & compartment == "medial_load"]
  pL <- as.numeric(table(factor(la$state, levels = retained))) / nrow(la)
  pM <- as.numeric(table(factor(me$state, levels = retained))) / nrow(me)
  muL <- vapply(retained, function(s) mean(la[state == s, PIEZO1]), numeric(1))
  muM <- vapply(retained, function(s) mean(me[state == s, PIEZO1]), numeric(1))
  if (any(!is.finite(c(muL, muM)))) stop("Retained state absent within donor/compartment")
  delta <- sum(pL * muL) - sum(pM * muM)
  within_L <- sum(pL * (muL - muM))
  composition_L <- sum(muM * (pL - pM))
  within_M <- sum(pM * (muL - muM))
  composition_M <- sum(muL * (pL - pM))
  within_symmetric <- (within_L + within_M) / 2
  composition_symmetric <- (composition_L + composition_M) / 2
  decomp_rows[[dd]] <- data.table(
    donor = dd,
    excluded_states = paste(low_recall, collapse = ";"),
    retained_states = paste(retained, collapse = ";"),
    n_lateral = nrow(la), n_medial = nrow(me), total_delta = delta,
    within_lateral_reference = within_L,
    composition_lateral_reference = composition_L,
    within_medial_reference = within_M,
    composition_medial_reference = composition_M,
    within_symmetric = within_symmetric,
    composition_symmetric = composition_symmetric,
    within_share_symmetric = within_symmetric / delta,
    composition_share_symmetric = composition_symmetric / delta,
    within_share_lateral_reference = within_L / delta,
    within_share_medial_reference = within_M / delta,
    decomposition_error = abs(within_symmetric + composition_symmetric - delta),
    medial_higher_within = within_symmetric < 0,
    both_components_contribute = sign(within_symmetric) == sign(composition_symmetric) &&
      sign(within_symmetric) == sign(delta)
  )
  direction_rows[[dd]] <- data.table(
    donor = dd, state = retained,
    n_lateral = as.numeric(table(factor(la$state, levels = retained))),
    n_medial = as.numeric(table(factor(me$state, levels = retained))),
    mean_lateral = muL, mean_medial = muM,
    delta_lateral_minus_medial = muL - muM
  )
}
decomp <- rbindlist(decomp_rows)
directions <- rbindlist(direction_rows)
stopifnot(max(decomp$decomposition_error) < 1e-10)
fwrite(decomp, file.path(OUT, "05_HighRecall_State_Sensitivity.csv"))
fwrite(directions, file.path(OUT, "05_HighRecall_State_Directions.csv"))
fwrite(recall[estimate < 0.60], file.path(OUT, "05_Excluded_LowRecall_States.csv"))

# R-only Figure S15. The hero panels expose the nested confusion matrix and the
# full-pipeline permutation null; validation panels show state recall, assignment
# margin, and the retained-state decomposition sensitivity.
pal <- c(EC = "#4477AA", FC = "#EE6677", HomC = "#228833", HTC = "#CCBB44",
         preHTC = "#AA3377", ProC = "#66CCEE", RegC = "#999999")
theme_set(theme_bw(base_size = 8.5) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", size = 9),
                  legend.title = element_text(size = 8)))

p_a <- ggplot(conf, aes(pred, truth, fill = row_fraction)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * row_fraction)), size = 2.25) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1),
                      name = "Row fraction") +
  labs(title = "a  Fully nested LODO confusion", x = "Predicted state", y = "Reference state") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_b <- ggplot(data.table(balanced_accuracy = null_ba), aes(balanced_accuracy)) +
  geom_histogram(bins = 35, fill = "#BDBDBD", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = observed$balanced_accuracy, color = "#B2182B", linewidth = 0.7) +
  annotate("text", x = observed$balanced_accuracy, y = Inf,
           label = sprintf("Observed BA = %.3f\nP = %.4f", observed$balanced_accuracy, perm_p),
           vjust = 1.2, hjust = 1.05, size = 2.5, color = "#B2182B") +
  labs(title = "b  Complete-pipeline permutation null", x = "Balanced accuracy", y = "Permutations")

p_c <- ggplot(recall, aes(truth, estimate, fill = truth)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0.60, linetype = 2, color = "#B2182B") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
  labs(title = "c  Held-out state recall", x = NULL, y = "Recall") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

margin_plot <- assignment[, .(margin = median(margin), q1 = quantile(margin, 0.25),
                              q3 = quantile(margin, 0.75)), by = state]
p_d <- ggplot(margin_plot, aes(state, margin, color = state)) +
  geom_linerange(aes(ymin = q1, ymax = q3), linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = pal, guide = "none") +
  labs(title = "d  Transfer-score separation", x = NULL, y = "Top-1 minus top-2 score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dec_long <- melt(decomp, id.vars = "donor",
                 measure.vars = c("within_symmetric", "composition_symmetric"),
                 variable.name = "component", value.name = "contribution")
dec_long[, donor := factor(donor, levels = sort(unique(donor)))]
dec_long[, component := factor(component,
  levels = c("within_symmetric", "composition_symmetric"),
  labels = c("Within-state", "Composition"))]
p_e <- ggplot(dec_long, aes(donor, contribution, fill = component)) +
  geom_col(position = "stack", width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  scale_fill_manual(values = c("Within-state" = "#2166AC", "Composition" = "#B2182B")) +
  labs(title = sprintf("e  High-recall-state sensitivity (excluded: %s)",
                       ifelse(length(low_recall), paste(low_recall, collapse = ", "), "none")),
       x = "Donor", y = "Lateral-minus-medial contribution", fill = NULL)

fig_s15 <- (p_a | p_b) / (p_c | p_d) / p_e +
  plot_annotation(title = "Reference-internal validation and sensitivity analyses for cross-dataset state transfer")

base <- file.path(FIG, "Figure_S15_Nested_StateTransfer_Validation")
ggsave(paste0(base, ".svg"), fig_s15, width = 183, height = 210, units = "mm",
       device = svglite::svglite)
ggsave(paste0(base, ".pdf"), fig_s15, width = 183, height = 210, units = "mm",
       device = cairo_pdf)
ggsave(paste0(base, ".png"), fig_s15, width = 183, height = 210, units = "mm",
       dpi = 300, device = ragg::agg_png)
ggsave(paste0(base, ".tiff"), fig_s15, width = 183, height = 210, units = "mm",
       dpi = 600, compression = "lzw", device = ragg::agg_tiff)

message(sprintf("Observed nested BA %.6f; macro-F1 %.6f; null95 %.6f; P %.6f",
                observed$balanced_accuracy, observed$macro_f1, null95, perm_p))
message(sprintf("Low-recall states: %s", ifelse(length(low_recall), paste(low_recall, collapse = ", "), "none")))
message("Nested validation and high-recall sensitivity outputs complete")
