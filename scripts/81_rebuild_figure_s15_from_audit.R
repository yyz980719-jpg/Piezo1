suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

source("scripts/_bootstrap.R")
OUT <- PIEZO1_METADATA
PRIOR <- PIEZO1_METADATA
FIG <- PIEZO1_FIGURES
qc <- fread(file.path(OUT, "04_Nested_LODO_QC.csv"))
conf <- fread(file.path(OUT, "04_LODO_confusion_matrix.csv"))
null <- fread(file.path(OUT, "04_Nested_LODO_PermutationNull.csv"))
decomp <- fread(file.path(OUT, "05_HighRecall_State_Sensitivity.csv"))
assignment <- fread(file.path(PRIOR, "01_GSE152805_StateAssignmentConfidence.csv"))
recall <- qc[metric == "recall", .(truth = state, estimate)]
states <- sort(unique(recall$truth))
ba <- qc[metric == "balanced_accuracy", estimate]
perm_p <- qc[metric == "balanced_accuracy", permutation_p]
low_recall <- recall[estimate < 0.60, truth]

pal <- c(EC = "#4477AA", FC = "#EE6677", HomC = "#228833", HTC = "#CCBB44",
         preHTC = "#AA3377", ProC = "#66CCEE", RegC = "#999999")
theme_set(theme_bw(base_size = 8.5) +
            theme(panel.grid.minor = element_blank(),
                  plot.title = element_text(face = "bold", size = 9),
                  legend.title = element_text(size = 8)))

p_a <- ggplot(conf, aes(pred, truth, fill = row_fraction)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * row_fraction)), size = 2.25) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1), name = "Row fraction") +
  labs(title = "a  Fully nested LODO confusion", x = "Predicted state", y = "Reference state") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_b <- ggplot(null, aes(balanced_accuracy)) +
  geom_histogram(bins = 35, fill = "#BDBDBD", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = ba, color = "#B2182B", linewidth = 0.7) +
  annotate("text", x = ba, y = Inf,
           label = sprintf("Observed BA = %.3f\nP = %.4f", ba, perm_p),
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
ggsave(paste0(base, ".svg"), fig_s15, width = 183, height = 210, units = "mm", device = svglite::svglite)
ggsave(paste0(base, ".pdf"), fig_s15, width = 183, height = 210, units = "mm", device = cairo_pdf)
ggsave(paste0(base, ".png"), fig_s15, width = 183, height = 210, units = "mm", dpi = 300, device = ragg::agg_png)
ggsave(paste0(base, ".tiff"), fig_s15, width = 183, height = 210, units = "mm", dpi = 600,
       compression = "lzw", device = ragg::agg_tiff)
