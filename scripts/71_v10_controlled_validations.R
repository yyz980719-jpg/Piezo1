# PIEZO1 OA V10 controlled validation workflow
# Implements SOP P0-01 to P0-04 without changing the frozen raw data.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
})

set.seed(20260831)
source("scripts/_bootstrap.R")
ROOT <- PIEZO1_ROOT
OUT <- PIEZO1_METADATA
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

theme_set(theme_bw(base_size = 8.5) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", size = 9),
                  strip.background = element_rect(fill = "grey94")))
pal <- c(EC="#4477AA", FC="#EE6677", HomC="#228833", HTC="#CCBB44",
         preHTC="#AA3377", ProC="#66CCEE", RegC="#BBBBBB")

softmax <- function(x) {
  x <- x - max(x)
  z <- exp(x)
  z / sum(z)
}

balanced_accuracy <- function(truth, pred, lev) {
  rr <- sapply(lev, function(s) {
    ix <- truth == s
    if (!any(ix)) return(NA_real_)
    mean(pred[ix] == s)
  })
  mean(rr, na.rm = TRUE)
}

macro_f1 <- function(truth, pred, lev) {
  ff <- sapply(lev, function(s) {
    tp <- sum(truth == s & pred == s)
    fp <- sum(truth != s & pred == s)
    fn <- sum(truth == s & pred != s)
    pr <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
    rc <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
    if ((pr + rc) == 0) 0 else 2 * pr * rc / (pr + rc)
  })
  mean(ff)
}

# -----------------------------------------------------------------------------
# Frozen reference marker construction (identical to scripts 33/39/52)
# -----------------------------------------------------------------------------
seu <- readRDS("data/scrna/GSE104782_seurat.rds")
md0 <- as.data.table(seu@meta.data, keep.rownames = "cell")
md0[, subtype := as.character(subtype)]
md0[, patient := as.character(patient)]
d0 <- tryCatch(GetAssayData(seu, layer = "data"),
               error = function(e) GetAssayData(seu, slot = "data"))
lin0 <- as.matrix(expm1(d0))
subtypes <- sort(unique(md0$subtype))
mean_full <- sapply(subtypes, function(s) rowMeans(lin0[, md0$subtype == s, drop = FALSE]))
marker_rank <- lapply(subtypes, function(s) {
  other <- rowMeans(mean_full[, setdiff(subtypes, s), drop = FALSE])
  fc <- (mean_full[, s] + 1e-6) / (other + 1e-6)
  fc[mean_full[, s] < 0.05] <- 0
  data.table(state = s, gene = names(sort(fc, decreasing = TRUE))[1:100],
             fold_change = sort(fc, decreasing = TRUE)[1:100])
})
marker_rank <- rbindlist(marker_rank)
sig_genes <- unique(marker_rank$gene)

# -----------------------------------------------------------------------------
# P0-01A: leave-one-donor-out state transfer in GSE104782
# -----------------------------------------------------------------------------
pred_parts <- list()
for (dd in sort(unique(md0$patient))) {
  train <- md0$patient != dd
  test <- md0$patient == dd
  cent <- sapply(subtypes, function(s) rowMeans(lin0[, train & md0$subtype == s, drop = FALSE]))
  sig <- cent[sig_genes, , drop = FALSE]
  sig_z <- t(scale(t(sig))); sig_z[is.na(sig_z)] <- 0
  q <- as.matrix(d0[sig_genes, test, drop = FALSE])
  qz <- t(scale(t(q))); qz[is.na(qz)] <- 0
  W <- t(qz) %*% sig_z
  pred <- subtypes[max.col(W, ties.method = "first")]
  pred_parts[[dd]] <- data.table(cell = md0$cell[test], donor = dd,
                                 truth = md0$subtype[test], pred = pred)
}
lod <- rbindlist(pred_parts)
conf <- as.data.table(table(factor(lod$truth, levels = subtypes),
                           factor(lod$pred, levels = subtypes)))
setnames(conf, c("truth", "pred", "N"))
rec <- lod[, .(recall = mean(pred == truth), n = .N), by = truth]
ba_obs <- balanced_accuracy(lod$truth, lod$pred, subtypes)
f1_obs <- macro_f1(lod$truth, lod$pred, subtypes)
nperm <- 1000L
null_ba <- numeric(nperm)
for (i in seq_len(nperm)) {
  ptruth <- lod[, sample(truth), by = donor]$V1
  null_ba[i] <- balanced_accuracy(ptruth, lod$pred, subtypes)
}
null95 <- unname(quantile(null_ba, 0.95, type = 8))
metrics <- rbind(
  data.table(metric = "balanced_accuracy", estimate = ba_obs, null_95th = null95,
             pass = ba_obs > null95),
  data.table(metric = "macro_F1", estimate = f1_obs, null_95th = NA_real_, pass = NA),
  rec[, .(metric = paste0("recall_", truth), estimate = recall,
          null_95th = NA_real_, pass = NA)]
)
fwrite(lod, file.path(OUT, "01_LODO_CellPredictions.csv"))
fwrite(conf, file.path(OUT, "01_LODO_ConfusionMatrix.csv"))
fwrite(data.table(iteration = seq_len(nperm), balanced_accuracy = null_ba),
       file.path(OUT, "01_LODO_PermutationNull.csv"))
fwrite(metrics, file.path(OUT, "01_StateTransfer_QC.csv"))

# -----------------------------------------------------------------------------
# Load GSE152805 and reproduce frozen projection scores
# -----------------------------------------------------------------------------
raw_dir <- "data/GSE152805_raw"
keep <- c("OA_oLT_113","OA_oLT_116","OA_oLT_118",
          "OA_MT_113","OA_MT_116","OA_MT_118")
counts <- NULL
meta <- data.table()
for (s in keep) {
  fm <- file.path(raw_dir, list.files(raw_dir, pattern = paste0("^GSM[0-9]+_", s, "\\.matrix")))
  fb <- file.path(raw_dir, list.files(raw_dir, pattern = paste0("^GSM[0-9]+_", s, "\\.barcodes")))
  fg <- file.path(raw_dir, list.files(raw_dir, pattern = paste0("^GSM[0-9]+_", s, "\\.genes")))
  M <- readMM(gzfile(fm))
  bc <- read.delim(gzfile(fb), header = FALSE, stringsAsFactors = FALSE)[, 1]
  gn <- read.delim(gzfile(fg), header = FALSE, stringsAsFactors = FALSE)
  dimnames(M) <- list(gn[, 2], paste(s, bc, sep = "_"))
  counts <- if (is.null(counts)) M else cbind(counts, M)
  meta <- rbind(meta, data.table(cell = colnames(M),
                                compartment = ifelse(grepl("oLT", s),
                                                     "lateral_nonload", "medial_load"),
                                donor = sub(".*_", "", s)))
}
nc <- Matrix::colSums(counts)
nf <- Matrix::colSums(counts > 0)
mtg <- grep("^MT-", rownames(counts), value = TRUE)
mito <- if (length(mtg)) Matrix::colSums(counts[mtg, , drop = FALSE]) / nc else rep(0, ncol(counts))
qc <- nf >= 200 & nc >= 500 & mito < 0.25
counts <- counts[, qc, drop = FALSE]
meta <- meta[qc]
if (any(duplicated(rownames(counts)))) counts <- counts[!duplicated(rownames(counts)), ]
lib <- Matrix::colSums(counts)
expr <- as.matrix(Matrix::t(log1p(Matrix::t(counts) / lib * 1e4)))

sig_full <- mean_full[sig_genes, , drop = FALSE]
sig_z_full <- t(scale(t(sig_full))); sig_z_full[is.na(sig_z_full)] <- 0
common <- intersect(sig_genes, rownames(expr))
zz <- t(scale(t(expr[common, , drop = FALSE]))); zz[is.na(zz)] <- 0
W <- t(zz) %*% sig_z_full[common, , drop = FALSE]
colnames(W) <- subtypes
ord <- t(apply(W, 1, order, decreasing = TRUE))
meta[, state := subtypes[ord[, 1]]]
meta[, max_weight := W[cbind(seq_len(nrow(W)), ord[, 1])]]
meta[, margin := W[cbind(seq_len(nrow(W)), ord[, 1])] - W[cbind(seq_len(nrow(W)), ord[, 2])]]
probs <- t(apply(W, 1, function(v) softmax(as.numeric(scale(v)))))
meta[, normalized_entropy := -rowSums(probs * log(probs + 1e-15)) / log(ncol(probs))]
meta[, PIEZO1 := if ("PIEZO1" %in% rownames(expr)) as.numeric(expr["PIEZO1", ]) else NA_real_]
fwrite(meta, file.path(OUT, "01_GSE152805_StateAssignmentConfidence.csv"))

# donor x predicted-state heatmap for up to top 10 reference markers/state
top_markers <- marker_rank[, head(.SD, 10), by = state]
heat_genes <- unique(intersect(top_markers$gene, rownames(expr)))
hm <- rbindlist(lapply(sort(unique(meta$donor)), function(dd) {
  rbindlist(lapply(subtypes, function(s) {
    ii <- which(meta$donor == dd & meta$state == s)
    if (!length(ii)) return(NULL)
    data.table(donor = dd, predicted_state = s, gene = heat_genes,
               mean_log_expression = rowMeans(expr[heat_genes, ii, drop = FALSE]))
  }))
}))
hm[, z := as.numeric(scale(mean_log_expression)), by = gene]
hm[is.na(z), z := 0]
fwrite(hm, file.path(OUT, "01_DonorPredictedState_MarkerHeatmap.csv"))

# -----------------------------------------------------------------------------
# P0-02 plus P0-01B: reference-sensitive decomposition and confidence trimming
# -----------------------------------------------------------------------------
decompose_one <- function(m) {
  out <- list()
  for (dd in sort(unique(m$donor))) {
    la <- m[donor == dd & compartment == "lateral_nonload"]
    me <- m[donor == dd & compartment == "medial_load"]
    pL <- as.numeric(table(factor(la$state, levels = subtypes))) / nrow(la)
    pM <- as.numeric(table(factor(me$state, levels = subtypes))) / nrow(me)
    muL <- sapply(subtypes, function(s) mean(la$PIEZO1[la$state == s]))
    muM <- sapply(subtypes, function(s) mean(me$PIEZO1[me$state == s]))
    if (anyNA(c(muL, muM))) stop("A state is absent within a donor/compartment after trimming")
    delta <- sum(pL * muL) - sum(pM * muM)
    WL <- sum(pL * (muL - muM)); CL <- sum(muM * (pL - pM))
    WM <- sum(pM * (muL - muM)); CM <- sum(muL * (pL - pM))
    WS <- (WL + WM) / 2; CS <- (CL + CM) / 2
    out[[dd]] <- data.table(donor = dd, total_delta = delta,
                            within_L_reference = WL, composition_L_reference = CL,
                            within_M_reference = WM, composition_M_reference = CM,
                            within_symmetric = WS, composition_symmetric = CS,
                            error_L = abs(WL + CL - delta),
                            error_M = abs(WM + CM - delta),
                            error_symmetric = abs(WS + CS - delta),
                            same_sign_symmetric = sign(WS) == sign(CS),
                            within_share_symmetric = ifelse(sign(WS) == sign(CS) && delta != 0,
                                                            WS / delta, NA_real_),
                            within_share_min = ifelse(sign(WL) == sign(CL) && delta != 0,
                                                      min(WL / delta, WM / delta), NA_real_),
                            within_share_max = ifelse(sign(WL) == sign(CL) && delta != 0,
                                                      max(WL / delta, WM / delta), NA_real_))
  }
  rbindlist(out)
}

trim_levels <- c(0, 0.10, 0.20)
dec_all <- list(); within_all <- list()
for (tr in trim_levels) {
  mm <- copy(meta)
  if (tr > 0) {
    mm[, cutoff := quantile(margin, probs = tr, type = 8), by = .(donor, compartment)]
    mm <- mm[margin > cutoff]
  }
  dec <- decompose_one(mm)
  dec[, trim_fraction := tr]
  dec_all[[as.character(tr)]] <- dec
  rows <- list()
  for (dd in sort(unique(mm$donor))) for (s in subtypes) {
    a <- mm[donor == dd & compartment == "lateral_nonload" & state == s, PIEZO1]
    b <- mm[donor == dd & compartment == "medial_load" & state == s, PIEZO1]
    rows[[length(rows) + 1]] <- data.table(trim_fraction = tr, donor = dd, state = s,
      n_lateral = length(a), n_medial = length(b),
      mean_lateral = mean(a), mean_medial = mean(b),
      delta_lateral_minus_medial = mean(a) - mean(b))
  }
  within_all[[as.character(tr)]] <- rbindlist(rows)
}
dec_all <- rbindlist(dec_all)
within_all <- rbindlist(within_all)
stopifnot(max(dec_all$error_L, dec_all$error_M, dec_all$error_symmetric) < 1e-10)
fwrite(dec_all, file.path(OUT, "02_Decomposition_Sensitivity.csv"))
fwrite(within_all, file.path(OUT, "01_LowConfidence_Sensitivity.csv"))

# -----------------------------------------------------------------------------
# P0-03: axis definition and donor-within quartiles
# -----------------------------------------------------------------------------
homeo_requested <- c("COL2A1","SOX9","ACAN","CHAD","HAPLN1")
hyp_requested <- c("COL10A1","RUNX2","IBSP","ALPL","MMP13","SPP1","POSTN")
homeo <- intersect(homeo_requested, rownames(d0))
hyp <- intersect(hyp_requested, rownames(d0))
score_set <- function(gs) as.numeric(scale(colMeans(d0[gs, , drop = FALSE])))
axis_score <- score_set(hyp) - score_set(homeo)
pz_ref <- as.numeric(d0["PIEZO1", ])
axis_cells <- copy(md0[, .(cell, patient, subtype)])
axis_cells[, axis_score := axis_score]
axis_cells[, PIEZO1 := pz_ref]
axis_cells[, quartile := {
  rr <- frank(axis_score, ties.method = "average") / .N
  fifelse(rr <= 0.25, "Q1", fifelse(rr > 0.75, "Q4", "middle"))
}, by = patient]
axis_donor <- axis_cells[quartile %chin% c("Q1", "Q4"),
                         .(n_cells = .N, mean_PIEZO1 = mean(PIEZO1),
                           mean_axis = mean(axis_score)), by = .(patient, quartile)]
axis_wide <- dcast(axis_donor, patient ~ quartile, value.var = "mean_PIEZO1")
axis_wide[, delta_Q4_minus_Q1 := Q4 - Q1]
tt_axis <- t.test(axis_wide$Q4, axis_wide$Q1, paired = TRUE)
wt_axis <- wilcox.test(axis_wide$Q4, axis_wide$Q1, paired = TRUE, exact = FALSE)
ct_axis <- suppressWarnings(cor.test(axis_cells$PIEZO1, axis_cells$axis_score, method = "spearman"))
axis_def <- rbind(
  data.table(item = "homeostatic_markers_requested", value = paste(homeo_requested, collapse = "; "),
             detail = "full prespecified set before intersection with the expression matrix"),
  data.table(item = "homeostatic_markers", value = paste(homeo, collapse = "; "),
             detail = "higher values decrease the axis score"),
  data.table(item = "homeostatic_alias_resolution", value = "CRTL1 -> HAPLN1",
             detail = "CRTL1 is a historical alias of HAPLN1 and is not an additional marker"),
  data.table(item = "hypertrophic_markers_requested", value = paste(hyp_requested, collapse = "; "),
             detail = "full prespecified set before intersection with the expression matrix"),
  data.table(item = "hypertrophic_markers", value = paste(hyp, collapse = "; "),
             detail = "higher values increase the axis score"),
  data.table(item = "hypertrophic_markers_absent", value = paste(setdiff(hyp_requested, hyp), collapse = "; "),
             detail = "requested genes absent from the normalized reference matrix"),
  data.table(item = "PIEZO1_in_marker_set", value = as.character("PIEZO1" %chin% c(homeo, hyp)),
             detail = "target-gene leave-out is unnecessary because PIEZO1 is absent"),
  data.table(item = "scaling", value = "cell-level gene-set mean; each set standardized across all reference cells",
             detail = "axis = z(mean hypertrophic markers) - z(mean homeostatic markers)"),
  data.table(item = "quartiles", value = "donor-within rank quartiles",
             detail = "Q1 <=25th percentile; Q4 >75th percentile within each donor"),
  data.table(item = "cell_level_spearman", value = sprintf("rho=%.6f; P=%.6g", ct_axis$estimate, ct_axis$p.value),
             detail = "descriptive cell-level association"),
  data.table(item = "donor_paired_Q4_minus_Q1", value = sprintf("mean=%.6f; paired t P=%.6g; Wilcoxon P=%.6g",
             mean(axis_wide$delta_Q4_minus_Q1), tt_axis$p.value, wt_axis$p.value),
             detail = sprintf("n=%d donors", nrow(axis_wide)))
)
fwrite(axis_def, file.path(OUT, "03_AxisDefinition.csv"))
fwrite(axis_cells, file.path(OUT, "03_AxisCells_DonorQuartiles.csv"))
fwrite(axis_wide, file.path(OUT, "03_AxisDonor_Q1Q4.csv"))

# -----------------------------------------------------------------------------
# P0-04: independent two-row manual dot-product check and units statement
# -----------------------------------------------------------------------------
manual_query <- c(1.25, -0.50)
manual_centroid <- c(0.80, 1.40)
manual_terms <- manual_query * manual_centroid
manual_score <- sum(manual_terms)
stopifnot(identical(as.numeric(crossprod(manual_query, manual_centroid)), manual_score))
unit_lines <- c(
  "Programme-score unit audit",
  "==========================",
  "Frozen implementation: W = t(query_gene_z) %*% reference_centroid_gene_z.",
  "Each cell-state score is therefore an unconstrained similarity-projection score.",
  "A donor-level compartment contrast is the difference between mean cell projection scores.",
  "The values are arbitrary but reproducible score units; they are neither proportions nor percentage points.",
  "",
  "Independent two-row manual dot-product test:",
  sprintf("query z = [%.2f, %.2f]", manual_query[1], manual_query[2]),
  sprintf("centroid z = [%.2f, %.2f]", manual_centroid[1], manual_centroid[2]),
  sprintf("terms = [%.2f, %.2f]", manual_terms[1], manual_terms[2]),
  sprintf("manual sum = %.2f; crossprod = %.2f; PASS", manual_score,
          as.numeric(crossprod(manual_query, manual_centroid)))
)
writeLines(unit_lines, file.path(OUT, "04_ProgrammeScore_Units.txt"))

# -----------------------------------------------------------------------------
# Figure S15: state-transfer validation and sensitivity
# -----------------------------------------------------------------------------
conf_plot <- copy(conf)
conf_plot[, row_total := sum(N), by = truth]
conf_plot[, prop := N / row_total]
pA <- ggplot(conf_plot, aes(pred, truth, fill = prop)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * prop)), size = 2.0) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1), name = "Row\nproportion") +
  labs(x = "Predicted state", y = "Published state", title = "LODO state-transfer confusion") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

null_dt <- data.table(x = null_ba)
pB <- ggplot(null_dt, aes(x)) +
  geom_histogram(bins = 35, fill = "grey78", colour = "white") +
  geom_vline(xintercept = null95, linetype = "dotted", colour = "#CC6677") +
  geom_vline(xintercept = ba_obs, linewidth = 0.8, colour = "#2166AC") +
  annotate("text", x = ba_obs, y = Inf, vjust = 1.3, hjust = 1.05,
           label = sprintf("Observed BA %.3f", ba_obs), size = 2.4, colour = "#2166AC") +
  labs(x = "Balanced accuracy", y = "Permutations", title = "Donor-preserving permutation null")

pC <- ggplot(meta, aes(state, margin, fill = state)) +
  geom_violin(scale = "width", colour = NA, alpha = 0.75) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.25) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(trans = "log1p") +
  labs(x = "Transferred state", y = "Top-1 minus top-2 score (log1p scale)", title = "Assignment margin") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

top3 <- unique(top_markers[, head(gene, 3), by = state]$V1)
hmplot <- hm[gene %chin% top3]
hmplot[, donor_state := paste0("D", donor, ".", predicted_state)]
pD <- ggplot(hmplot, aes(donor_state, gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Gene-wise z") +
  labs(x = "Donor.predicted state", y = NULL, title = "Reference-marker expression") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5.5),
        axis.text.y = element_text(size = 5.5))

decplot <- melt(dec_all[, .(donor, trim_fraction, within_symmetric, composition_symmetric)],
                id.vars = c("donor", "trim_fraction"), variable.name = "component")
decplot[, component := factor(component,
  levels = c("within_symmetric", "composition_symmetric"),
  labels = c("within-state expression", "composition"))]
pE <- ggplot(decplot, aes(factor(trim_fraction * 100), value, fill = component)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  facet_wrap(~ donor, nrow = 1) +
  scale_fill_manual(values = c("#4477AA", "#CC6677"), name = NULL) +
  labs(x = "Lowest-margin cells removed within donor x compartment (%)",
       y = "Contribution to lateral - medial PIEZO1",
       title = "Symmetric decomposition after confidence trimming") +
  theme(legend.position = "top")

design <- "
AB
CD
EE
"
fig <- pA + pB + pC + pD + pE +
  plot_layout(design = design, heights = c(1.00, 1.25, 0.80)) +
  plot_annotation(tag_levels = "a",
                  title = "Supplementary Figure S15 | Validation of cross-dataset state transfer") &
  theme(plot.tag = element_text(face = "bold", size = 9))

base <- file.path(OUT, "Figure_S15_StateTransfer_Validation")
ggsave(paste0(base, ".svg"), fig, width = 183, height = 220, units = "mm", bg = "white")
ggsave(paste0(base, ".pdf"), fig, width = 183, height = 220, units = "mm", bg = "white",
       device = cairo_pdf)
ggsave(paste0(base, ".tiff"), fig, width = 183, height = 220, units = "mm", dpi = 600,
       compression = "lzw", bg = "white")
ggsave(paste0(base, ".png"), fig, width = 183, height = 220, units = "mm", dpi = 400,
       bg = "white")

# concise run summary used by the document-edit branch logic
summary_lines <- c(
  sprintf("LODO balanced accuracy: %.6f", ba_obs),
  sprintf("Permutation null 95th percentile: %.6f", null95),
  sprintf("LODO gate: %s", ifelse(ba_obs > null95, "PASS", "FAIL")),
  sprintf("LODO macro-F1: %.6f", f1_obs),
  sprintf("Maximum decomposition identity error: %.3g", max(dec_all$error_L, dec_all$error_M, dec_all$error_symmetric)),
  sprintf("Symmetric within-state components at 0%% trim: %s",
          paste(sprintf("D%s=%+.6f", dec_all[trim_fraction == 0]$donor,
                        dec_all[trim_fraction == 0]$within_symmetric), collapse = "; ")),
  sprintf("Symmetric within-state shares at 0%% trim: %s",
          paste(sprintf("D%s=%s", dec_all[trim_fraction == 0]$donor,
                        ifelse(is.na(dec_all[trim_fraction == 0]$within_share_symmetric), "not applicable",
                               sprintf("%.1f%%", 100 * dec_all[trim_fraction == 0]$within_share_symmetric))),
                collapse = "; ")),
  sprintf("Donor-within axis Q4-Q1 PIEZO1 mean difference: %.6f", mean(axis_wide$delta_Q4_minus_Q1)),
  sprintf("Donor-level paired t P: %.6g; Wilcoxon P: %.6g", tt_axis$p.value, wt_axis$p.value)
)
writeLines(summary_lines, file.path(OUT, "00_ControlledValidation_Summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
