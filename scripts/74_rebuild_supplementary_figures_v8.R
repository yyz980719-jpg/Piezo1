suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(svglite); library(ragg); library(scales)
})

source("scripts/_bootstrap.R")
ROOT <- PIEZO1_ROOT
OUT <- PIEZO1_FIGURES
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

pal <- c(blue="#247B9A", blue2="#8BB5C7", orange="#C9563B", gold="#E2A44F",
         grey="#999999", light="#E7E7E7", green="#4D9F94", dark="#2B2B2B")
theme_set(theme_classic(base_size = 8, base_family = "Arial") +
            theme(plot.title = element_text(face = "bold", size = 9),
                  plot.tag = element_text(face = "bold", size = 10),
                  strip.background = element_blank(), strip.text = element_text(face = "bold")))

save_pub <- function(p, stem, width_mm = 183, height_mm = 105) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite(file.path(OUT, paste0(stem, ".svg")), width = w, height = h); print(p); dev.off()
  cairo_pdf(file.path(OUT, paste0(stem, ".pdf")), width = w, height = h, family = "Arial"); print(p); dev.off()
  agg_tiff(file.path(OUT, paste0(stem, ".tiff")), width = w, height = h,
           units = "in", res = 600, compression = "lzw"); print(p); dev.off()
  agg_png(file.path(OUT, paste0(stem, ".png")), width = w, height = h,
          units = "in", res = 400); print(p); dev.off()
}

# -----------------------------------------------------------------------------
# Supplementary Figure S3: remove target-gene self-correlation and rerank.
# -----------------------------------------------------------------------------
lat <- fread("results/07_GO_LateralUp.csv")
med <- fread("results/07_GO_MedialUp.csv")
terms_lat <- c("regulation of lymphocyte activation", "fatty acid metabolic process", "response to insulin")
terms_med <- c("mitotic nuclear division", "skeletal system development", "extracellular matrix organization")
go <- rbind(
  lat[Description %chin% terms_lat, .SD[1], by = Description][, direction := "Lateral"],
  med[Description %chin% terms_med, .SD[1], by = Description][, direction := "Medial"]
)
go[, score := fifelse(direction == "Lateral", 1, -1) * -log10(p.adjust)]
go[, Description := factor(Description, levels = Description[order(score)])]
p3a <- ggplot(go, aes(score, Description)) +
  geom_vline(xintercept = 0, colour = pal["grey"], linewidth = .35) +
  geom_segment(aes(x = 0, xend = score, yend = Description, colour = direction), linewidth = .65) +
  geom_point(aes(size = Count, fill = direction), shape = 21, colour = "white") +
  scale_colour_manual(values = c(Lateral = unname(pal["blue"]), Medial = unname(pal["orange"])), guide = "none") +
  scale_fill_manual(values = c(Lateral = unname(pal["blue"]), Medial = unname(pal["orange"])), guide = "none") +
  labs(x = "Signed -log10(FDR)", y = NULL, size = "Genes",
       title = "Opposing regional programmes",
       subtitle = "Blue: lateral enriched; orange: medial enriched") +
  theme(legend.position = "bottom")

corr <- fread("results/07_PIEZO1_delta_correlation.csv")
setnames(corr, names(corr)[1:2], c("gene", "rho"))
corr <- corr[gene != "PIEZO1"][order(-rho)][1:14]
corr[, gene := factor(gene, levels = rev(gene))]
p3b <- ggplot(corr, aes(rho, gene)) +
  geom_col(fill = pal["green"], width = .72) +
  labs(x = "Spearman rho with PIEZO1 donor delta", y = NULL,
       title = "Top exploratory co-variation genes",
       subtitle = "PIEZO1 excluded from its own correlation ranking") +
  coord_cartesian(xlim = c(0, max(corr$rho) * 1.05))

figS3 <- (p3a | p3b) + plot_layout(widths = c(1.08, 1)) + plot_annotation(tag_levels = "a")
save_pub(figS3, "Figure_S3_RegionalEnrichment_Covariation_v8", 183, 88)
fwrite(corr, metadata_path("SourceData_FigureS3b_target_excluded.csv"))

# -----------------------------------------------------------------------------
# Supplementary Figure S7: academically neutral title and explicit score unit.
# -----------------------------------------------------------------------------
pd <- fread("results/39_gse152805_perdonor.csv")
setnames(pd, names(pd)[1], "measure")
prog <- pd[measure != "PIEZO1"]
pz <- pd[measure == "PIEZO1"]
subtypes <- c("EC", "FC", "HomC", "HTC", "preHTC", "ProC", "RegC")
dl <- melt(prog, id.vars = "measure", measure.vars = subtypes,
           variable.name = "programme", value.name = "delta")
setnames(dl, "measure", "donor")
dl[, donor := factor(donor, levels = c("113", "116", "118"))]
dl[, programme := factor(programme, levels = subtypes)]
pzlong <- melt(pz, id.vars = "measure", measure.vars = subtypes,
               variable.name = "donor", value.name = "delta")
pzlong <- data.table(donor = factor(subtypes[1:3], levels = c("113", "116", "118")), delta = NA_real_)
pzlong <- data.table(donor = factor(c("113", "116", "118"), levels = c("113", "116", "118")),
                     delta = as.numeric(unlist(pz[1, ..subtypes]))[1:3])

p7a <- ggplot(dl, aes(programme, donor, fill = delta)) +
  geom_tile(colour = "white", linewidth = .45) +
  scale_fill_gradient2(low = pal["orange"], mid = "white", high = pal["blue"], midpoint = 0,
                       name = "Lateral - medial") +
  labs(x = "Reference programme", y = "Donor", title = "All programme directions",
       subtitle = "Unconstrained similarity-projection score units") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
p7b <- ggplot(pzlong, aes(delta, donor)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_segment(aes(x = 0, xend = delta, yend = donor), linewidth = 1.05, colour = pal["orange"]) +
  geom_point(size = 3.1, colour = pal["orange"]) +
  labs(x = "PIEZO1 difference (lateral - medial)", y = NULL,
       title = "PIEZO1: medial higher in all three donors")
figS7 <- (p7a | p7b) + plot_layout(widths = c(1, 1.05)) +
  plot_annotation(tag_levels = "a",
                  title = "Supplementary Figure S7 | GSE152805 donor-level directional concordance")
save_pub(figS7, "Figure_S7_GSE152805_DirectionalConcordance_v8", 183, 80)
fwrite(dl, metadata_path("SourceData_FigureS7_programme_projection_scores.csv"))

# -----------------------------------------------------------------------------
# Supplementary Figure S12: symmetric decomposition, no pooled percentage.
# -----------------------------------------------------------------------------
within <- fread(metadata_path("01_LowConfidence_Sensitivity.csv"))[trim_fraction == 0]
within[, state := factor(state, levels = c("EC", "HomC", "HTC", "RegC", "FC", "preHTC", "ProC"))]
within[, donor := factor(donor, levels = c("113", "116", "118"))]
p12a <- ggplot(within, aes(state, delta_lateral_minus_medial, colour = donor, size = pmin(n_lateral, n_medial))) +
  geom_hline(yintercept = 0, linewidth = .35) +
  geom_point(position = position_dodge(.58), alpha = .9) +
  scale_colour_manual(values = c("113" = unname(pal["blue"]), "116" = unname(pal["orange"]), "118" = unname(pal["gold"]))) +
  scale_size_area(max_size = 5) +
  labs(x = NULL, y = "Lateral - medial PIEZO1", colour = "Donor",
       size = "Smaller cell count", title = "Within-state expression differences") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

comp <- fread("results/52_state_composition.csv")
comp[, state := factor(state, levels = levels(within$state))]
comp[, donor := factor(donor, levels = c("113", "116", "118"))]
p12b <- ggplot(comp, aes(state, delta_prop, fill = donor)) +
  geom_hline(yintercept = 0, linewidth = .35) + geom_col(position = "dodge") +
  scale_fill_manual(values = c("113" = unname(pal["blue"]), "116" = unname(pal["orange"]), "118" = unname(pal["gold"]))) +
  labs(x = NULL, y = "Lateral - medial state proportion", fill = "Donor",
       title = "State-composition differences") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

dec <- fread(metadata_path("02_Decomposition_Sensitivity.csv"))[trim_fraction == 0]
dm <- melt(dec[, .(donor, within_symmetric, composition_symmetric)], id.vars = "donor",
           variable.name = "component", value.name = "value")
dm[, component := factor(component, levels = c("within_symmetric", "composition_symmetric"),
                         labels = c("Within-state expression", "Composition"))]
p12c <- ggplot(dm, aes(factor(donor), value, fill = component)) +
  geom_hline(yintercept = 0, linewidth = .35) + geom_col() +
  geom_point(data = dec, aes(factor(donor), total_delta), inherit.aes = FALSE,
             shape = 18, size = 2.6) +
  scale_fill_manual(values = c("Within-state expression" = unname(pal["blue"]), "Composition" = unname(pal["light"]))) +
  labs(x = "Donor", y = "Contribution to lateral - medial PIEZO1", fill = NULL,
       title = "Symmetric decomposition") + theme(legend.position = "top")

figS12 <- p12a / p12b / p12c + plot_layout(heights = c(1, 1, .92)) +
  plot_annotation(tag_levels = "a",
    title = "Supplementary Figure S12 | Within-state expression and compositional contributions")
save_pub(figS12, "Figure_S12_SymmetricDecomposition_v8", 183, 190)

cat("Supplementary figures S3, S7 and S12 exported to", OUT, "\n")
