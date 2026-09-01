#!/usr/bin/env Rscript
# 19_celltype_enrichment.R
# Phase 2-C: GWAS cell-type enrichment (lightweight MAGMA/LDSC alternative)
#   OA GWAS (FinnGen R12 M13 arthrosis, full) gene-level p via nearest_genes
#   x scRNA marker genes (GSE324993, chondrocyte subclusters)
#   Hypergeometric enrichment of OA-significant genes within each cell-type marker set.
suppressMessages({
  library(data.table)
})
source("scripts/_bootstrap.R")
ROOT <- PIEZO1_ROOT
GW   <- data_path("mr", "finngen_R12_M13_ARTHROSIS.gz")
MK   <- file.path(ROOT, "data/scrna/GSE324993_marker_genes.txt.gz")
OUTC <- file.path(ROOT, "results/19_celltype_enrichment.csv")
OUTF <- file.path(ROOT, "figs/fig68_celltype_enrichment.png")

cat("Reading OA GWAS (gene-level via nearest_genes)...\n"); flush.console()
gw <- fread(GW, select = c("nearest_genes", "pval"), showProgress = FALSE)
gw <- gw[!is.na(pval) & nearest_genes != "" & !grepl("^intergenic", nearest_genes, ignore.case = TRUE)]
# expand comma-separated nearest genes
genes <- strsplit(gw$nearest_genes, ",", fixed = TRUE)
gdt   <- data.table(gene = unlist(genes), pval = rep(gw$pval, lengths(genes)))
gmin  <- gdt[, .(p = min(pval)), by = gene]
N     <- nrow(gmin)                       # background = genes with GWAS data
cat(sprintf("  background genes N=%d\n", N)); flush.console()

thr <- 1e-5
O   <- gmin[p < thr, unique(gene)]
q   <- length(O)
cat(sprintf("  OA-significant genes (p<%.0e): q=%d\n", thr, q)); flush.console()

mk <- fread(MK)
mk_sig <- mk[p_val_adj < 0.05]
clusters <- unique(mk_sig$Cluster)
cat(sprintf("  marker clusters: %s\n", paste(clusters, collapse = ", "))); flush.console()

rows <- list()
for (cl in clusters) {
  S <- unique(mk_sig[Cluster == cl, Genes])
  ov <- intersect(S, O)
  k <- length(S); m <- length(ov)
  # hypergeometric: P(X >= m)
  p_h <- phyper(m - 1, q, N - q, k, lower.tail = FALSE)
  exp <- k * q / N
  OR  <- if (m > 0 && (k - m) > 0 && (q - m) > 0 && (N - q - k + m) > 0)
           (m / (k - m)) / ((q - m) / (N - q - k + m)) else NA_real_
  rows[[cl]] <- data.table(cluster = cl, n_markers = k, n_OA_sig = m,
                           expected = round(exp, 2), OR = round(OR, 2),
                           p_hyperg = p_h, p_fdr = NA_real_)
}
out <- rbindlist(rows)
out[, p_fdr := p.adjust(p_hyperg, method = "BH")]
setorder(out, -n_OA_sig)
fwrite(out, OUTC)
cat("Wrote", OUTC, "\n"); print(out)

# ---- figure: -log10(p) per cluster, sized by n_OA_sig ----
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressMessages(library(ggplot2))
  p <- ggplot(out, aes(x = reorder(cluster, -p_hyperg), y = -log10(p_hyperg),
                      size = n_OA_sig, fill = -log10(p_hyperg))) +
    geom_point(shape = 21, color = "black") +
    scale_fill_gradient(low = "steelblue", high = "red") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = "OA GWAS gene enrichment by chondrocyte subcluster (Phase 2-C)",
         subtitle = sprintf("background N=%d genes; OA-sig(p<1e-5) q=%d; hypergeometric", N, q),
         x = "scRNA cluster (marker gene set)", y = "-log10(hypergeometric p)") +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(OUTF, p, width = 8, height = 5, dpi = 150)
  cat("Wrote", OUTF, "\n")
} else {
  cat("ggplot2 not available; CSV only\n")
}
cat("DONE Phase2-C\n")
