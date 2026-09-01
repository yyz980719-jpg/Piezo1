# ==========================================================================
# Nature-style main figures for the PIEZO1 x OA manuscript.
# Backend: R (ggplot2 + patchwork). Data are read from stored result tables.
# FIGURE 4  Cross-dataset replication of the chondrocyte-programme contrast
#           and the divergent direction of PIEZO1 transcription
# ==========================================================================

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

source("scripts/_bootstrap.R")
RES  <- PIEZO1_RESULTS
WORK <- PIEZO1_METADATA
OUT  <- PIEZO1_FIGURES
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Nature typography / theme -------------------------------------------
theme_nature <- function(base = 6.5) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(
      axis.line  = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text  = element_text(size = base - 0.5, colour = "black"),
      axis.title = element_text(size = base, colour = "black"),
      legend.title    = element_text(size = base - 0.3),
      legend.text     = element_text(size = base - 0.7),
      legend.key.size = unit(3, "mm"),
      legend.background = element_blank(),
      strip.text = element_text(size = base - 0.3, face = "bold", hjust = 0),
      strip.background = element_blank(),
      plot.title = element_text(size = base + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base - 0.6, colour = COL_NEUTRAL, hjust = 0),
      plot.tag = element_text(size = base + 2.5, face = "bold"),
      plot.tag.location = "plot",
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 1, 1)
    )
}
COL_NEUTRAL <- "#595959"
COL_LIGHT   <- "#C9C9C9"
COL_SIGNAL  <- "#1F6F8B"   # programme contrast / primary
COL_ACCENT  <- "#C1553B"   # PIEZO1 / diverging quantity
COL_WARM    <- "#E0A458"
COL_POS     <- "#1F6F8B"
COL_NEG     <- "#C1553B"
theme_set(theme_nature())

save_nature <- function(p, name, w_mm = 183, h_mm = 120, dpi = 600) {
  w <- w_mm / 25.4; h <- h_mm / 25.4
  svglite::svglite(file.path(OUT, paste0(name, ".svg")), width = w, height = h)
  print(p); dev.off()
  grDevices::cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h,
                       family = "Arial"); print(p); dev.off()
  ragg::agg_tiff(file.path(OUT, paste0(name, ".tiff")), width = w, height = h,
                 units = "in", res = dpi, compression = "lzw"); print(p); dev.off()
  ragg::agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h,
                units = "in", res = 400); print(p); dev.off()
  cat("[saved]", name, "\n")
}

# ==========================================================================
# data
# ==========================================================================
cs   <- fread(file.path(RES, "51_compositional_sensitivity_full.csv"))
pd   <- fread(file.path(RES, "39_gse152805_perdonor.csv"))
setnames(pd, 1, "donor")
pd   <- pd[donor != "PIEZO1"]
wst  <- fread(file.path(RES, "52_within_state_per_donor.csv"))[available == TRUE]
dec  <- fread(file.path(RES, "52_kitagawa_decomposition.csv"))
ovr  <- fread(file.path(RES, "52_overall_direction.csv"))
axd  <- fread(file.path(RES, "40_P0_donorlevel_axis_perdonor.csv"))
rep_g <- fread(file.path(RES, "35_replication_genes.csv"))

FAMMAP <- c(L0 = "Raw", L1 = "CLR", L2 = "ALR", L3 = "Pairwise", L4 = "GS score")
cs[, family := FAMMAP[substr(layer, 1, 2)]]

# ---- build contrast / isolated-HomC table --------------------------------
bld <- rbindlist(lapply(split(cs, list(cs$dataset, cs$layer)), function(d) {
  ds <- d$dataset[1]; ly <- d$layer[1]; is_alr <- grepl("^L2", ly)
  h <- d[programme == "HomC"]; p <- d[programme == "preHTC"]
  dir1 <- d[programme %in% c("log(HomC/preHTC)", "HomC_minus_preHTC")]
  r <- list()
  # contrast: ALR is reference-invariant; a missing coordinate is identically zero
  if (nrow(h) && nrow(p))
    r[[length(r) + 1]] <- data.table(dataset = ds, layer = ly, quantity = "contrast",
                                     delta = h$delta - p$delta)
  else if (is_alr && !nrow(h) && nrow(p))
    r[[length(r) + 1]] <- data.table(dataset = ds, layer = ly, quantity = "contrast",
                                     delta = -p$delta)
  else if (is_alr && nrow(h) && !nrow(p))
    r[[length(r) + 1]] <- data.table(dataset = ds, layer = ly, quantity = "contrast",
                                     delta = h$delta)
  if (nrow(dir1))
    r[[length(r) + 1]] <- data.table(dataset = ds, layer = ly, quantity = "contrast",
                                     delta = dir1$delta)
  if (nrow(h))
    r[[length(r) + 1]] <- data.table(dataset = ds, layer = ly, quantity = "HomC alone",
                                     delta = h$delta)
  rbindlist(r)
}))
bld <- unique(bld, by = c("dataset", "layer", "quantity", "delta"))
bld[, family  := FAMMAP[substr(layer, 1, 2)]]
bld[, dataset := factor(dataset,
                        levels = c("GSE51588(LT-MT)", "GSE57218(Pres-OA)"),
                        labels = c("GSE51588 bone (LT - MT)", "GSE57218 cartilage (preserved - OA)"))]
bld[, quantity := factor(quantity, levels = c("contrast", "HomC alone"),
                         labels = c("HomC minus preHTC contrast", "HomC alone"))]
bld[, family := factor(family, levels = c("Raw", "CLR", "ALR", "Pairwise", "GS score"))]

stopifnot(bld[quantity == "HomC minus preHTC contrast", all(delta > 0)])
cat("contrast positive in",
    bld[quantity == "HomC minus preHTC contrast", .N], "of",
    bld[quantity == "HomC minus preHTC contrast", .N], "layer instances per dataset\n")
cat("HomC alone:", bld[quantity == "HomC alone", sum(delta > 0)], "positive /",
    bld[quantity == "HomC alone", sum(delta <= 0)], "non-positive\n")

# ==========================================================================
# (a) evidence hierarchy  -- schematic
# ==========================================================================
tier <- data.table(
  y0 = c(74, 48, 48, 20), y1 = c(96, 70, 70, 42),
  x0 = c(30,  2, 52, 30), x1 = c(98, 48, 98, 98),
  lab  = c("GSE104782", "GSE57218", "GSE152805", "GSE51588"),
  sub  = c("single-cell OA cartilage\n1,464 cells, 10 donors",
           "paired preserved / OA cartilage\n33 donors",
           "paired lateral / medial cartilage\n3 donors",
           "paired medial / lateral subchondral bone\n20 OA + 5 non-OA donors"),
  role = c("Programme definition", "Paired replication",
           "Paired replication", "Regional discovery"),
  tier_col = c(COL_WARM, COL_SIGNAL, COL_SIGNAL, COL_NEUTRAL)
)
role_box <- data.table(y0 = c(74, 48, 20), y1 = c(96, 70, 42), x0 = 0, x1 = 27,
                       txt = c("Defines the\nseven states",
                               "Replication of the\nprogramme contrast",
                               "Regional discovery;\ncross-tissue check"))
arrows <- data.table(x = 62, xend = 62, y = c(74, 48), yend = c(70, 42))

pa <- ggplot() +
  geom_rect(data = tier, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = "#F4F4F2", colour = "black", linewidth = 0.3) +
  geom_rect(data = tier, aes(xmin = x0, xmax = x0 + 1.6, ymin = y0, ymax = y1),
            fill = tier$tier_col, colour = NA) +
  geom_rect(data = role_box, aes(xmin = x0, xmax = x1, ymin = y0, ymax = y1),
            fill = NA, colour = COL_LIGHT, linewidth = 0.3, linetype = "dashed") +
  geom_text(data = role_box, aes(x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = txt),
            size = 1.75, colour = COL_NEUTRAL, lineheight = 0.85) +
  geom_text(data = tier, aes(x = x0 + 4, y = y1 - 5, label = lab),
            hjust = 0, size = 2.35, fontface = "bold", colour = "black") +
  geom_text(data = tier, aes(x = x0 + 4, y = y1 - 13, label = role),
            hjust = 0, size = 1.85, colour = COL_SIGNAL) +
  geom_text(data = tier, aes(x = x0 + 4, y = y0 + 6.5, label = sub),
            hjust = 0, vjust = 0, size = 1.7, colour = COL_NEUTRAL, lineheight = 0.9) +
  geom_segment(data = arrows, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
               linewidth = 0.35, colour = "black") +
  scale_x_continuous(limits = c(-1, 100)) +
  scale_y_continuous(limits = c(16, 99)) +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), plot.margin = margin(1, 1, 1, 1))

# ==========================================================================
# (b) five analytic layers
# ==========================================================================
fam_lab <- bld[, .(n = uniqueN(layer)), by = .(family = as.character(family))]
fam_lab[, family := factor(family, levels = levels(bld$family))]
setorder(fam_lab, family)
fam_lab <- bld[, .(n = uniqueN(layer)), by = family][order(family)]
fam_lab[, lbl := paste0(as.character(family), "\n(", n, ifelse(n > 1, " rules)", " rule)"))]

pos_frac <- bld[, .(n = .N, npos = sum(delta > 0)), by = .(quantity, family)]
pb <- ggplot(bld, aes(x = family, y = delta, colour = dataset)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.35) +
  geom_point(position = position_dodge(width = 0.62), size = 1.15, alpha = 0.9,
             shape = 16) +
  scale_colour_manual(values = c("GSE51588 bone (LT - MT)" = COL_SIGNAL,
                                 "GSE57218 cartilage (preserved - OA)" = COL_ACCENT),
                      name = NULL) +
  scale_x_discrete(labels = setNames(fam_lab$lbl, as.character(fam_lab$family))) +
  facet_wrap(~ quantity, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "Estimate (layer- and dataset-specific units)") +
  theme(axis.text.x = element_text(size = 5.6, colour = "black"),
        legend.position = "top", legend.margin = margin(0, 0, 1, 0),
        legend.text = element_text(size = 5.8),
        strip.text = element_text(size = 6.2, face = "bold"))

# ==========================================================================
# (c) programme direction per GSE152805 donor
# ==========================================================================
PROG <- c("EC", "HomC", "HTC", "RegC", "FC", "preHTC", "ProC")
pd_m <- melt(pd[donor %in% c("113", "116", "118")],
             id.vars = "donor", variable.name = "programme", value.name = "delta")
pd_m[, donor := factor(donor, levels = c("113", "116", "118"))]
pd_m[, programme := factor(programme, levels = rev(PROG))]
pd_m[, dir := ifelse(delta > 0, "lateral higher", "medial higher")]
agree <- pd_m[, .(npos = sum(delta > 0)), by = programme]
agree[, lbl := ifelse(npos == 3, "3/3", ifelse(npos == 0, "3/3", "2/3"))]

pc <- ggplot(pd_m, aes(x = donor, y = programme, fill = dir)) +
  geom_tile(colour = "white", linewidth = 0.4, width = 0.86, height = 0.86) +
  geom_text(aes(label = ifelse(delta > 0, "+", "\u2212")),
            size = 3.4, colour = "white", fontface = "bold") +
  geom_text(data = agree, aes(x = 3.85, y = programme, label = lbl),
            inherit.aes = FALSE, size = 1.75, colour = "black", hjust = 0) +
  annotate("text", x = 3.85, y = 7.75, label = "concordant",
           size = 1.7, colour = COL_NEUTRAL, hjust = 0, fontface = "italic") +
  scale_fill_manual(values = c("lateral higher" = COL_POS, "medial higher" = COL_NEG),
                    name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.35))) +
  labs(x = "Donor", y = NULL) +
  theme(legend.position = "top", legend.margin = margin(0, 0, 1, 0),
        legend.text = element_text(size = 5.6),
        axis.text = element_text(colour = "black"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())

# ==========================================================================
# (d) within-state PIEZO1 difference by state and donor
# ==========================================================================
wst[, donor := factor(donor, levels = c("113", "116", "118"))]
wst[, state := factor(state, levels = c("EC", "HomC", "HTC", "RegC", "FC", "preHTC", "ProC"))]
n_neg <- wst[, .(n = .N, nneg = sum(delta_lateral_minus_medial < 0)), by = donor]

pd_d <- ggplot(wst, aes(x = state, y = delta_lateral_minus_medial,
                        colour = donor, shape = donor)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.35) +
  geom_point(position = position_dodge(width = 0.6), size = 1.35, alpha = 0.95) +
  scale_colour_manual(values = c("113" = COL_SIGNAL, "116" = COL_ACCENT, "118" = COL_WARM),
                      name = "Donor") +
  scale_shape_manual(values = c(16, 17, 15), name = "Donor") +
  scale_y_continuous(labels = label_number(style_negative = "minus")) +
  labs(x = NULL, y = "Lateral minus medial\nPIEZO1 (within state)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        legend.position = "top", legend.margin = margin(0, 0, 1, 0),
        legend.text = element_text(size = 5.6))

# ==========================================================================
# (e) Kitagawa decomposition
# ==========================================================================
dec_long <- rbindlist(list(
  dec[, .(donor = as.character(donor), component = "Within state", value = within_state_component)],
  dec[, .(donor = as.character(donor), component = "Composition", value = composition_component)]
))
dec_long[, component := factor(component, levels = c("Within state", "Composition"))]
dec_lab <- copy(dec)[, `:=`(donor = as.character(donor),
                            ylab = total_delta_lateral_minus_medial + 0.004,
                            txt = sprintf("%.0f%% within", pct_within))]
agg <- dec_long[, .(value = sum(value)), by = component]
agg_within <- -sum(dec$within_state_component) / -sum(dec$total_delta_lateral_minus_medial) * 100

pe <- ggplot(dec_long, aes(x = donor, y = value, fill = component)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.35) +
  geom_col(width = 0.6, colour = NA) +
  geom_point(data = dec[, .(donor = as.character(donor),
                            value = total_delta_lateral_minus_medial)],
             inherit.aes = FALSE, aes(x = donor, y = value),
             shape = 18, size = 2.1, colour = "black") +
  geom_text(data = dec_lab, aes(x = donor, y = ylab, label = txt),
            inherit.aes = FALSE, size = 1.75, colour = "black", hjust = 0.3) +
  scale_fill_manual(values = c("Within state" = COL_SIGNAL, "Composition" = COL_LIGHT),
                    name = NULL) +
  scale_y_continuous(labels = label_number(style_negative = "minus")) +
  labs(x = "Donor", y = "Lateral minus medial PIEZO1\n(decomposed)") +
  theme(legend.position = "top", legend.margin = margin(0, 0, 1, 0),
        legend.text = element_text(size = 5.6),
        axis.text = element_text(colour = "black")) +
  annotate("text", x = 0.45, y = -0.008,
           label = sprintf("pooled: %.0f%% within-state", agg_within),
           size = 1.7, colour = COL_NEUTRAL, hjust = 0, fontface = "italic")

# ==========================================================================
# (f) PIEZO1 direction on dataset-specific scales
#     convention: lateral / preserved minus medial / OA-affected
# ==========================================================================
rho  <- axd$rho
rho_m <- mean(rho); rho_se <- sd(rho) / sqrt(length(rho))
rho_ci <- rho_m + c(-1, 1) * qt(0.975, length(rho) - 1) * rho_se

pz5 <- -(-0.5035281644)                       # lateral - medial, OA donors
pz5_ci <- c(0.260878589061453, 0.746177739738548)
pz7 <- rep_g[gene == "PIEZO1", delta_preserved_minus_OA]
pz7_t <- rep_g[gene == "PIEZO1", t]
pz7_se <- abs(pz7 / pz7_t); pz7_ci <- pz7 + c(-1, 1) * qt(0.975, 32) * pz7_se
pz8 <- ovr$overall_delta_lateral_minus_medial

f_dat <- rbindlist(list(
  data.table(src = "GSE104782", est = rho_m, lo = rho_ci[1], hi = rho_ci[2],
             unit = "\u03c1 with hypertrophic-minus-homeostatic axis", n = "10 donors",
             tissue = "Cartilage"),
  data.table(src = "GSE152805", est = mean(pz8), lo = NA, hi = NA,
             unit = "lateral minus medial mean difference", n = "3 donors",
             tissue = "Cartilage"),
  data.table(src = "GSE57218", est = pz7, lo = pz7_ci[1], hi = pz7_ci[2],
             unit = "preserved minus OA-affected difference", n = "33 donors",
             tissue = "Cartilage"),
  data.table(src = "GSE51588", est = pz5, lo = pz5_ci[1], hi = pz5_ci[2],
             unit = "lateral minus medial difference, OA donors", n = "20 donors",
             tissue = "Subchondral bone")
))
f_dat[, ylab := paste0(src, "  \u00b7  ", unit, "\n", n)]
f_dat[, ylab := factor(ylab, levels = rev(ylab))]
f_pts <- data.table(y = pz8, ylab = f_dat[src == "GSE152805", ylab])

pf <- ggplot() +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.35) +
  geom_errorbar(data = f_dat[!is.na(lo)], aes(xmin = lo, xmax = hi, y = ylab),
                width = 0.10, linewidth = 0.4, colour = "black") +
  geom_point(data = f_dat[!is.na(lo)], aes(x = est, y = ylab, fill = est > 0),
             shape = 21, size = 2.2, colour = "black", stroke = 0.25) +
  geom_point(data = f_pts, aes(x = y, y = ylab, fill = y > 0),
             shape = 21, size = 1.7, colour = "black", stroke = 0.25,
             position = position_nudge(y = c(-0.14, 0, 0.14))) +
  geom_text(data = f_dat[!is.na(lo)], aes(x = est, y = ylab,
            label = sprintf("%.3f", est)), vjust = -1.9, size = 1.7,
            colour = "black") +
  scale_fill_manual(values = c("TRUE" = COL_POS, "FALSE" = COL_NEG), guide = "none") +
  labs(x = "Estimate on the dataset's own scale (95% CI where estimable)", y = NULL) +
  theme(axis.text.y = element_text(colour = "black", size = 5.9, hjust = 0),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 6))

# ==========================================================================
# assemble
# ==========================================================================
fig4 <- pa / pb / (pc | pd_d | pe) / pf +
  plot_layout(heights = c(0.85, 1.55, 1.15, 1.0)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))

save_nature(fig4, "Fig4_programme_contrast_PIEZO1_direction", w_mm = 183, h_mm = 205)

# ---- source data ----------------------------------------------------------
fwrite(bld, file.path(OUT, "sourcedata_Fig4b_layers.csv"))
fwrite(pd_m, file.path(OUT, "sourcedata_Fig4c_programme_direction.csv"))
fwrite(wst, file.path(OUT, "sourcedata_Fig4d_within_state.csv"))
fwrite(dec, file.path(OUT, "sourcedata_Fig4e_kitagawa.csv"))
fwrite(f_dat, file.path(OUT, "sourcedata_Fig4f_piezo1_direction.csv"))
cat("\nFig 4 numbers\n")
cat(sprintf("  rho  = %.3f [%.3f, %.3f]\n", rho_m, rho_ci[1], rho_ci[2]))
cat(sprintf("  GSE152805 lateral-medial: %s (mean %.3f)\n",
            paste(sprintf("%.3f", pz8), collapse = ", "), mean(pz8)))
cat(sprintf("  GSE57218 preserved-OA = %.3f [%.3f, %.3f]\n", pz7, pz7_ci[1], pz7_ci[2]))
cat(sprintf("  GSE51588 lateral-medial = %.3f [%.3f, %.3f]\n", pz5, pz5_ci[1], pz5_ci[2]))
cat(sprintf("  pooled Kitagawa within-state = %.1f%%\n", agg_within))
cat(sprintf("  within-state donors negative: %s\n",
            paste(sprintf("%s:%d/7", n_neg$donor, n_neg$nneg), collapse = " ")))

# ---- docx-sized render ----
D1 <- file.path(PIEZO1_METADATA, "figure4_docx_preview")
dir.create(D1, showWarnings = FALSE, recursive = TRUE)
ragg::agg_png(file.path(D1, "main_04.png"), width = 6.4, height = 6.4 * 205 / 183,
              units = "in", res = 300); print(fig4); dev.off()
cat("[docx] main_04.png\n")
