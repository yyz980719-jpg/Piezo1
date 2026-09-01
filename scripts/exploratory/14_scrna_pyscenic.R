# =====================================================================
# 14_scrna_pyscenic.R  —  Single-cell regulon activity (AUCell / SCENIC-lite)
# Upstream of PIEZO1 at single-cell resolution, using AUC-based regulon
# activity (the standard pySCENIC readout) instead of wmean.
# Input : data/scrna/GSE104782_seurat.rds (already processed, 7 subtypes)
# Method: dorothea_hs TF->target regulons; AUCell AUC per cell; then
#         correlate regulon AUC with PIEZO1 expression (per-cell + per-subtype)
# Output: figs/fig57_*.png, results/14_*.csv
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(dorothea); library(AUCell); library(data.table)
  library(ggplot2); library(cowplot); library(gridExtra)
})

source("scripts/_bootstrap.R")
so   <- readRDS("data/scrna/GSE104782_seurat.rds")
emat <- as.matrix(GetAssayData(so, layer = "data"))
gene_symbols <- rownames(emat)

# PIEZO1 per-cell expression
pz <- as.numeric(emat["PIEZO1", ])

# ---- build regulons from DoRothEA (tf -> target genes present in scRNA) ----
reg <- dorothea_hs
reg <- reg[!is.na(reg$mor) & reg$mor != 0, ]
reg <- reg[reg$target %in% gene_symbols, ]
reg <- reg[reg$tf %in% gene_symbols, ]
set.seed(1)
reg <- reg[sample(nrow(reg), min(nrow(reg), 200000)), ]   # speed
reglist <- split(reg$target, reg$tf)
reglist <- reglist[sapply(reglist, function(g) length(g) >= 15)]
cat(sprintf("[regulons] %d TFs with >=15 target genes in scRNA\n", length(reglist)))

# ---- AUCell rankings + AUC ----
rankings <- AUCell::AUCell_buildRankings(emat, plotStats = FALSE, verbose = FALSE)
regulonAUC <- AUCell::AUCell_calcAUC(reglist, rankings, aucMaxRank = 0.1*nrow(emat), verbose = FALSE)
auc <- AUCell::getAUC(regulonAUC)        # regulons (TF) x cells
auc <- auc[, colnames(emat)]

# ---- correlate each regulon AUC with PIEZO1 (per-cell) ----
sub_lab <- so$sub_lab
rho <- sapply(rownames(auc), function(tf)
  suppressWarnings(cor(auc[tf, ], pz, method = "spearman")))
pval <- sapply(rownames(auc), function(tf)
  suppressWarnings(cor.test(auc[tf, ], pz, method = "spearman")$p.value))
tab <- data.frame(tf = rownames(auc), rho = as.numeric(rho), p = as.numeric(pval))
tab$fdr <- p.adjust(tab$p, "BH")
tab <- tab[order(-abs(tab$rho)), ]
fwrite(tab, "results/14_scrna_regulon_PIEZO1_corr.csv")
cat(sprintf("[corr] regulons tested=%d | |rho|>=0.3 & FDR<0.05 : %d\n",
            nrow(tab), sum(abs(tab$rho)>=0.3 & tab$fdr<0.05)))

# cross-check with the 8 consensus TFs from bulk+single-cell target enrichment
consensus <- c("TCF7L1","GATAD2A","PRDM16","BCL6","ZNF92","ZNF853","PKNOX2","ATF6")
tab$consensus <- tab$tf %in% consensus

# ============ FIG 57: bar of |rho| (top 25), consensus highlighted ============
top <- tab[order(-abs(tab$rho)), ][1:25, ]
top$tf <- factor(top$tf, levels = top$tf[order(abs(top$rho))])
p57 <- ggplot(top, aes(x = tf, y = abs(rho), fill = consensus)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = c("FALSE"="#2166AC","TRUE"="#B2182B"),
                    name = "consensus TF\n(bulk+sc valid.)") +
  coord_flip() + theme_minimal(base_size = 10) +
  labs(title = "Fig S11A. Single-cell regulon activity correlated with PIEZO1",
       subtitle = "AUCell AUC of DoRothEA regulons vs PIEZO1 expression (Spearman, all cells)",
       x = NULL, y = "|Spearman rho|") +
  theme(legend.position = "bottom")
ggsave("figs/fig57_scrna_regulon_PIEZO1_bar.png", p57, width = 7, height = 7, dpi = 300)

# ============ FIG 58: UMAP of top regulon + PIEZO1 ============
top_tf <- tab$tf[which.max(abs(tab$rho))]
um <- as.data.frame(Embeddings(so, "umap"))
colnames(um) <- c("UMAP_1", "UMAP_2")
um$PIEZO1 <- pz
um$regulon <- as.numeric(auc[top_tf, ])
um$sub <- sub_lab
p58a <- ggplot(um, aes(UMAP_1, UMAP_2, color = PIEZO1)) + geom_point(size = 0.6) +
  scale_color_viridis_c(option = "magma") +
  ggtitle(paste0("PIEZO1 expression")) + theme_minimal(base_size = 9) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
p58b <- ggplot(um, aes(UMAP_1, UMAP_2, color = regulon)) + geom_point(size = 0.6) +
  scale_color_viridis_c(option = "viridis") +
  ggtitle(paste0(top_tf, " regulon AUC")) + theme_minimal(base_size = 9) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
p58 <- plot_grid(p58a, p58b, ncol = 2)
ggsave("figs/fig58_scrna_topregulon_umap.png", p58, width = 8, height = 4, dpi = 300)
cat(sprintf("[fig58] top regulon = %s (rho=%.3f)\n", top_tf, tab$rho[which.max(abs(tab$rho))]))

# ============ FIG 59: per-subtype rho for the 8 consensus TFs ============
cons_rows <- tab[tab$tf %in% consensus, ]
sub_types <- sort(unique(sub_lab))
sub_rho <- matrix(NA, nrow = length(consensus), ncol = length(sub_types),
                  dimnames = list(consensus, sub_types))
for (tf in consensus) {
  if (!tf %in% rownames(auc)) next
  for (s in sub_types) {
    idx <- sub_lab == s
    if (sum(idx) > 5)
      sub_rho[tf, s] <- suppressWarnings(cor(auc[tf, idx], pz[idx], method = "spearman"))
  }
}
sub_rho <- as.data.frame(sub_rho)
sub_rho$tf <- rownames(sub_rho)
melt <- data.table::melt(as.data.table(sub_rho), id.vars = "tf", variable.name = "subtype", value.name = "rho")
melt$subtype <- as.character(melt$subtype)
melt$tf <- factor(melt$tf, levels = consensus)
p59 <- ggplot(melt, aes(x = subtype, y = tf, fill = rho)) + geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       na.value = "grey90", name = "Spearman rho") +
  theme_minimal(base_size = 10) +
  labs(title = "Fig S11B. Regulon–PIEZO1 correlation by chondrocyte subtype",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figs/fig59_scrna_consensus_subtype_heatmap.png", p59, width = 8, height = 5, dpi = 300)

cat("\n[14] Done. figs fig57-59 + results/14_*.csv written.\n")
