suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

source("scripts/_bootstrap.R")

# Figure contract
# Core conclusion: the HomC-minus-preHTC programme contrast is directionally
# robust under all final analytic specifications, whereas PIEZO1 direction differs
# across state, compartment and tissue contexts. Magnitudes are not compared
# across transformations or datasets.
# Evidence chain: a, hierarchy; b, robustness matrix; c, three-donor
# programme direction; d, within-state PIEZO1; e, decomposition; f, direction-only
# cross-context synthesis.
# Archetype: quantitative grid. Backend: R only.
# Export: 183 x 205 mm, editable SVG/PDF, 600-dpi TIFF, 400-dpi PNG preview.

source(
  "scripts/fig_nature_4_sanitized.R",
  encoding = "UTF-8"
)

# The parent plotting script emits the L4 HomC-minus-preHTC contrast twice:
# once as HomC - preHTC and once from the explicit precomputed contrast row.
# They are algebraically identical but differ at floating-point precision.
# Keep one authoritative layer-level record, yielding 13 instances per dataset.
bld <- bld[, .SD[1L], by = .(dataset, layer, quantity)]

# Rebuild panel a as a compact evidence-tier strip. The parent panel was sized
# for a standalone export and its dense annotations overlap in the final grid.
tier_a <- data.table(
  x = 1:4,
  dataset = c("GSE104782", "GSE57218", "GSE152805", "GSE51588"),
  role = c(
    "Programme definition\n10 donors",
    "Paired-cartilage\nconcordance; 33 donors",
    "External projection support\n3 donors",
    "Regional discovery;\ncross-tissue similarity\n20 OA + 5 non-OA donors"
  ),
  group = factor(c("Definition", "Concordance", "Support", "Discovery"),
                 levels = c("Definition", "Concordance", "Support", "Discovery"))
)
pa <- ggplot(tier_a, aes(x = x, y = 1, fill = group)) +
  geom_tile(width = 0.94, height = 0.72, colour = "black", linewidth = 0.35) +
  geom_text(aes(label = dataset), y = 1.11, size = 2.7, fontface = "bold", colour = "white") +
  geom_text(aes(label = role), y = 0.89, size = 1.85, lineheight = 0.90, colour = "white") +
  scale_fill_manual(values = c(
    "Definition" = COL_WARM,
    "Concordance" = COL_SIGNAL,
    "Support" = "#3A8CA5",
    "Discovery" = COL_NEUTRAL
  ), guide = "none") +
  scale_x_continuous(NULL, breaks = NULL, limits = c(0.5, 4.5), expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, limits = c(0.55, 1.45), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1))

# Remove the clipped free-floating header from panel c; the adjacent 3/3 or
# 2/3 labels are retained and are defined in the legend as donor concordance.
pc <- ggplot(pd_m, aes(x = donor, y = programme, fill = dir)) +
  geom_tile(colour = "white", linewidth = 0.4, width = 0.86, height = 0.86) +
  geom_text(aes(label = ifelse(delta > 0, "+", "\u2212")),
            size = 3.4, colour = "white", fontface = "bold") +
  geom_text(data = agree, aes(x = 3.85, y = programme, label = lbl),
            inherit.aes = FALSE, size = 1.75, colour = "black", hjust = 0) +
  scale_fill_manual(values = c("lateral higher" = COL_POS, "medial higher" = COL_NEG),
                    name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.5, 1.35))) +
  labs(x = "Donor", y = NULL) +
  theme_nature(7.2) +
  theme(
    legend.position = "top", legend.margin = margin(0, 0, 1, 0),
    legend.text = element_text(size = 5.6),
    axis.text = element_text(colour = "black"),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.margin = margin(2, 5, 1, 1)
  )

DEST <- PIEZO1_FIGURES
dir.create(DEST, recursive = TRUE, showWarnings = FALSE)

# Panel b: conclusions by analytic layer, never cross-scale magnitudes.
matrix_b <- bld[, {
  vals <- delta
  if (all(vals > 0)) {
    status <- "Consistent +"
  } else if (all(vals < 0)) {
    status <- "Consistent -"
  } else {
    status <- "Reference-dependent"
  }
  list(status = status, label = sprintf("%d/%d +", sum(vals > 0), .N))
}, by = .(dataset, quantity, family)]

full_b <- CJ(
  dataset = levels(bld$dataset),
  quantity = levels(bld$quantity),
  family = levels(bld$family),
  unique = TRUE
)
matrix_b <- merge(full_b, matrix_b, by = c("dataset", "quantity", "family"), all.x = TRUE)
matrix_b[is.na(status), `:=`(status = "Not applicable", label = "N/A")]
matrix_b[, dataset := factor(dataset, levels = levels(bld$dataset))]
matrix_b[, quantity := factor(quantity, levels = rev(levels(bld$quantity)))]
matrix_b[, family := factor(family, levels = levels(bld$family))]

n_contrast <- bld[quantity == "HomC minus preHTC contrast", .N]
stopifnot(n_contrast == 26L)

pb_new <- ggplot(matrix_b, aes(family, quantity, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.55, width = 0.94, height = 0.88) +
  geom_text(aes(label = label), size = 2.05, fontface = "bold", colour = "white") +
  facet_wrap(~dataset, ncol = 2) +
  scale_fill_manual(values = c(
    "Consistent +" = COL_SIGNAL,
    "Consistent -" = COL_ACCENT,
    "Reference-dependent" = COL_WARM,
    "Not applicable" = COL_LIGHT
  ), name = NULL) +
  labs(
    x = "Analytic layer", y = NULL,
    subtitle = paste(strwrap(
      "The contrast remained positive under every final transformation, zero-handling rule and reference specification in both datasets.",
      width = 92), collapse = "\n")
  ) +
  theme_nature(7.2) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 6.2),
    axis.text.x = element_text(size = 6.2, colour = "black"),
    axis.text.y = element_text(size = 6.2, colour = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 6.7, face = "bold")
  )

# Remove the pooled percentage because reference-sensitive donor shares cross
# 50%; retain donor-specific values and the exact additive decomposition only.
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
  theme_nature(7.2) +
  theme(legend.position = "top", legend.margin = margin(0, 0, 1, 0),
        legend.text = element_text(size = 5.6),
        axis.text = element_text(colour = "black"))

# Panel f: categorical direction-only synthesis. Numerical estimates remain in
# the source-data table and legend; no common effect axis is drawn.
dir_f <- data.table(
  dataset = factor(
    c("GSE104782", "GSE152805", "GSE57218", "GSE51588"),
    levels = rev(c("GSE104782", "GSE152805", "GSE57218", "GSE51588"))
  ),
  context = c(
    "Across-state association",
    "Cartilage compartment",
    "Cartilage region",
    "Subchondral-bone compartment"
  ),
  direction = factor(c("Negative", "Negative", "Negative", "Positive"),
                     levels = c("Negative", "Positive")),
  statement = c(
    "PIEZO1 decreases toward hypertrophic pole; 9/10 donor correlations negative",
    "Lateral < medial PIEZO1; 3/3 donors",
    "Preserved < OA-affected PIEZO1; 33 paired donors",
    "Lateral > medial PIEZO1; 20 OA paired donors"
  )
)
dir_f[, row_label := paste0(as.character(dataset), "  |  ", context)]

pf_new <- ggplot(dir_f, aes(x = 1, y = dataset, fill = direction)) +
  geom_tile(width = 1.92, height = 0.76, colour = "white", linewidth = 0.7) +
  geom_text(aes(label = statement), colour = "white", fontface = "bold",
            size = 2.45, lineheight = 0.94) +
  scale_fill_manual(values = c("Negative" = COL_ACCENT, "Positive" = COL_SIGNAL), name = NULL) +
  scale_x_continuous(NULL, breaks = NULL, limits = c(0, 2)) +
  scale_y_discrete(labels = setNames(dir_f$row_label, as.character(dir_f$dataset))) +
  labs(y = NULL, subtitle = "Direction only; effect magnitudes are not compared across contexts") +
  theme_nature(7.2) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 6.2),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 6.2, colour = "black", face = "bold"),
    panel.grid = element_blank()
  )

fig4_final <- pa / pb_new / (pc | pd_d | pe) / pf_new +
  plot_layout(heights = c(0.88, 1.35, 1.08, 1.00)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9.5, face = "bold"))

save_final <- function(p, stem, w_mm = 183, h_mm = 205, dpi = 600) {
  w <- w_mm / 25.4
  h <- h_mm / 25.4
  svglite(file.path(DEST, paste0(stem, ".svg")), width = w, height = h)
  print(p)
  dev.off()
  cairo_pdf(file.path(DEST, paste0(stem, ".pdf")), width = w, height = h, family = "Arial")
  print(p)
  dev.off()
  agg_tiff(file.path(DEST, paste0(stem, ".tiff")), width = w, height = h,
           units = "in", res = dpi, compression = "lzw")
  print(p)
  dev.off()
  agg_png(file.path(DEST, paste0(stem, ".png")), width = w, height = h,
          units = "in", res = 400)
  print(p)
  dev.off()
}

stem <- "Figure4_programme_robustness_direction_v11"
save_final(fig4_final, stem)
fwrite(matrix_b, file.path(DEST, "SourceData_Figure4b_direction_robustness.csv"))
fwrite(dir_f[, .(dataset, context, direction, statement)],
       file.path(DEST, "SourceData_Figure4f_direction_only.csv"))
fwrite(bld[quantity == "HomC minus preHTC contrast"],
       file.path(DEST, "SourceData_Figure4_26_instances.csv"))

qa <- c(
  "Core conclusion: HomC-minus-preHTC is positive under every final specification in both datasets; PIEZO1 direction is context dependent.",
  "Backend: R only.",
  "Panel b: categorical direction/robustness matrix; no cross-transformation magnitude axis.",
  "Panel f: categorical direction-only rows; no cross-dataset magnitude axis.",
  "GSE152805: three-donor directional support only.",
  "Pooled percentage removed because donor-specific reference-sensitive shares cross 50%.",
  "Exports: SVG, PDF, 600-dpi TIFF, 400-dpi PNG; 183 x 205 mm."
)
writeLines(qa, file.path(DEST, "Figure4_QA_notes.txt"))
cat("Final Figure 4 written to", DEST, "\n")
