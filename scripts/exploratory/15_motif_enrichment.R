# =====================================================================
# 15_motif_enrichment.R  —  Sequence-level motif enrichment in the
# PIEZO1+ chondrocyte program (orthogonal evidence of upstream TFs)
# Input : data/scrna/GSE104782_seurat.rds  (PIEZO1+ program genes)
# Method: promoters (TSS +/-500bp, hg38) -> JASPAR2024 motif scanning
#         (motifmatchr) -> Fisher over-representation of each motif in
#         the PIEZO1+ program vs background (all expressed genes)
# Output: figs/fig60_*.png, results/15_motif_enrichment.csv
# =====================================================================
suppressPackageStartupMessages({
  library(Seurat); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg38); library(motifmatchr)
  library(JASPAR2024);   library(org.Hs.eg.db); library(data.table)
  library(SummarizedExperiment); library(ggplot2)
})

source("scripts/_bootstrap.R")
so   <- readRDS("data/scrna/GSE104782_seurat.rds")
emat <- as.matrix(GetAssayData(so, layer = "data"))
pz   <- as.numeric(emat["PIEZO1", ])

# ---- PIEZO1+ program: genes whose scRNA expression correlates with PIEZO1 ----
rho <- apply(emat, 1, function(g) suppressWarnings(cor(g, pz, method = "spearman")))
rho <- rho[!is.na(rho)]
expressed <- names(rho)[rho > 0 & !is.na(rho)]
# background = genes expressed in >=10% of cells
frac_expr <- rowMeans(emat > 0)
bg_sym <- names(frac_expr[frac_expr > 0.10])
bg_sym <- intersect(bg_sym, rownames(emat))
# program = top 250 positively correlated genes (excluding PIEZO1 itself)
prog <- names(sort(rho[rho > 0 & names(rho) != "PIEZO1"], decreasing = TRUE))[1:250]
prog <- intersect(prog, bg_sym)
cat(sprintf("[program] n=%d  [background] n=%d\n", length(prog), length(bg_sym)))

# ---- map SYMBOL -> ENTREZ for promoter extraction ----
sym2eg <- mapIds(org.Hs.eg.db, bg_sym, "ENTREZID", "SYMBOL")
eg2sym <- setNames(bg_sym, sym2eg[bg_sym])
eg_universe <- eg2sym[!is.na(eg2sym)]
universe_sym <- eg2sym[!is.na(eg2sym)]

# ---- promoters (TSS +/-500bp), standard chromosomes only ----
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gn   <- genes(txdb, filter = list(gene_id = names(universe_sym)))
seqlevelsStyle(gn) <- "UCSC"
gn   <- keepStandardChromosomes(gn, pruning.mode = "coarse")
gn   <- gn[!grepl("_", seqnames(gn))]            # drop alt/patched contigs
gn   <- gn[names(gn) %in% names(universe_sym)]   # keep mapped
pro  <- promoters(gn, upstream = 500, downstream = 100)
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, pro)
names(seqs) <- universe_sym[names(gn)]           # named by SYMBOL

# ---- motif scanning (MotifDb: human TF motifs, static annotation, no download) ----
suppressMessages({ library(MotifDb); library(TFBSTools); library(motifmatchr) })
mm  <- MotifDb
hs  <- mm[which(values(mm)$organism == "Hsapiens")]
hs  <- hs[!is.na(values(hs)$geneSymbol)]         # keep TF motifs with a symbol
make_pwm <- function(m, id, nm) {
  p <- m + 1e-4; p <- p / colSums(p); bg <- rep(0.25, 4); pwm <- log2(p / bg)
  PWMatrix(ID = id, name = nm, matrixClass = "PWM", profileMatrix = pwm)
}
lst <- list()
for (i in seq_along(hs))
  lst[[i]] <- tryCatch(make_pwm(hs[[i]], names(hs)[i], values(hs)$geneSymbol[i]), error = function(e) NULL)
lst <- Filter(Negate(is.null), lst)
pwml <- do.call(PWMatrixList, lst)
cat(sprintf("[motifs] scanning %d promoters vs %d human TF motifs (MotifDb)...\n",
            length(seqs), length(pwml)))
mres <- matchMotifs(pwml, seqs, out = "matches")
mat  <- tryCatch(assays(mres)[["motifMatches"]], error = function(e) NULL)
if (is.null(mat)) mat <- assay(mres)
mat  <- as.matrix(mat)                           # genes x motifs (counts)
# matchMotifs drops dimnames; restore from source order (subject = seqs, pwms = hs)
colnames(mat) <- names(hs)
rownames(mat) <- names(seqs)
pres <- mat > 0                                 # presence matrix
motif_names <- colnames(pres)
cat(sprintf("[motifs] matched motif columns = %d (regions=%d)\n", ncol(mat), nrow(mat)))
cat(sprintf("[DBG] motif_names len=%d | class(pres)=%s | nrow(pres)=%d\n",
            length(motif_names), class(pres), nrow(pres)))

# ---- over-representation (Fisher) for each motif ----
prog_set <- prog
rn <- rownames(pres)
idx_p <- rn %in% prog_set
idx_b <- rn %in% setdiff(names(seqs), prog_set)
res <- lapply(seq_along(motif_names), function(i) {
  m <- motif_names[i]
  in_p <- sum(pres[idx_p, m])
  in_b <- sum(pres[idx_b, m])
  ft <- fisher.test(matrix(c(in_p, sum(idx_p)-in_p, in_b, sum(idx_b)-in_b), 2, 2))
  data.frame(motif = m, n_target = in_p, n_bg = in_b,
             OR = ft$estimate, p = ft$p.value)
})
enr <- rbindlist(res)
cat(sprintf("[DBG] res len=%d | enr rows=%d | prog_set len=%d\n",
            length(res), nrow(enr), length(prog_set)))
enr$fdr <- p.adjust(enr$p, "BH")
# TF name from MotifDb motif id (column header is the MotifDb ID; symbol in rownames of values)
motif_sym <- setNames(as.character(values(hs)$geneSymbol), names(hs))
enr$tf  <- motif_sym[enr$motif]
enr <- enr[order(enr$fdr, -enr$OR), ]
fwrite(enr, "results/15_motif_enrichment.csv")
cat(sprintf("[motif] motifs tested=%d | FDR<0.05 : %d\n",
            nrow(enr), sum(enr$fdr < 0.05)))

# ---- FIG 60: top enriched motifs ----
top <- enr[enr$fdr < 0.05 & enr$n_target >= 3]
if (nrow(top) == 0) top <- enr[order(enr$fdr), ][1:min(25, nrow(enr)), ]
if (nrow(top) > 25) top <- top[1:25, ]
if (nrow(top) > 0) {
  top$lab <- paste0(top$tf, " (", top$motif, ")")
  top$lab <- factor(top$lab, levels = top$lab[order(-log10(top$fdr))])
  p60 <- ggplot(top, aes(x = lab, y = -log10(fdr), fill = n_target)) +
    geom_col(width = 0.8, color = "white") +
    scale_fill_gradient(low = "#DEEBF7", high = "#08519C", name = "n program genes") +
    coord_flip() + theme_minimal(base_size = 10) +
    labs(title = "Fig S12. TF motifs enriched in the PIEZO1+ chondrocyte program",
         subtitle = "MotifDb human TF motifs over-representation (Fisher, promoters TSS+/-500bp, hg38)",
         x = NULL, y = "-log10(FDR)") +
    theme(plot.title = element_text(face = "bold"))
  ggsave("figs/fig60_motif_enrichment.png", p60, width = 8, height = max(5, 0.32*nrow(top)+2), dpi = 300)
} else {
  cat("[fig60] no enriched motifs to plot\n")
}

cat("\n[15] Done. fig60 + results/15_motif_enrichment.csv written.\n")
