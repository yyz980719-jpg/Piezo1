# ============================================================
# 09_gsea.R — GSEA 秩次法通路富集分析 (PIEZO1 × OA)
# 四组基因排序 × MSigDB Hallmark (fgsea):
#   R1: GSE51588 软骨下骨 区域配对 (MT负重 vs LT非负重, 供体阻断)
#   R2: GSE51588 疾病主效应 (OA vs Normal, 校正区域)
#   R3: GSE55235 滑膜 OA vs Healthy (limma t)
#   R4: GSE82107 滑膜 OA vs Healthy (limma t)
#   R5: PIEZO1 共表达排序 (GSE51588 全 OA 样本 Spearman)
# ============================================================
source("scripts/_bootstrap.R")
suppressMessages({
  library(data.table); library(ggplot2); library(limma)
  library(GEOquery); library(fgsea); library(msigdbr)
  library(hgu133a.db); library(hgu133plus2.db)
  library(ggrepel)
})
figd <- "figs"; resd <- "results"
set.seed(42)

## ---------- Hallmark 基因集 ----------
m <- msigdbr(species = "Homo sapiens")
h <- m[grepl("^HALLMARK_", m$gs_name), ]
paths <- lapply(split(h$gene_symbol, h$gs_name), unique)
cat("Hallmark 通路数:", length(paths), "\n")

## ---------- 通用: GSEA + 气泡图 ----------
run_fgsea <- function(stats, label) {
  stats <- sort(stats[!is.na(stats)], decreasing = TRUE)
  fg <- fgsea(pathways = paths, stats = stats, minSize = 10, maxSize = 500,
              eps = 0, scoreType = "std")
  fg <- as.data.table(fg)[order(pval)]
  fg[, comparison := label]
  cat(sprintf("[%s] 显著通路 (FDR<0.05): %d / %d\n", label,
              sum(fg$padj < 0.05), nrow(fg)))
  fg
}

bubble <- function(fg, title, topn = 20, subtitle = NULL) {
  d <- as.data.table(fg)
  d <- d[padj < 0.25][order(-abs(NES))]
  if (nrow(d) > topn) d <- rbind(d[NES > 0][order(padj)][1:ceiling(topn/2)],
                                 d[NES < 0][order(padj)][1:ceiling(topn/2)])
  d <- unique(d, by = "pathway")
  d[, laby := gsub("HALLMARK_", "", pathway)]
  d[, laby := gsub("_", " ", laby)]
  d[, laby := factor(laby, levels = laby[order(NES)])]
  ggplot(d, aes(NES, laby, color = padj, size = size)) +
    geom_point(alpha = 0.9) +
    scale_color_gradientn(colors = c("#B2182B", "#EF8A62", "#FDDBC7", "#D1E5F0", "#67A9CF", "#2166AC"),
                           trans = "reverse", name = "FDR q") +
    scale_size_continuous(range = c(2, 7), name = "Gene count") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    labs(x = "Normalized enrichment score NES (right = pathway activation)", y = NULL,
         title = title, subtitle = subtitle) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          axis.text.y = element_text(size = 8.5))
}

## ============================================================
## ---------- GSE51588 软骨下骨 (Agilent GPL13497) ----------
## ============================================================
ann <- fread("data/GPL13497_annot_full.tsv", header = FALSE,
             col.names = c("probe", "symbol"))
ann <- ann[symbol != "" & !is.na(symbol)]
g5 <- getGEO(filename = "data/GSE51588_series_matrix.txt.gz", getGPL = FALSE)
e5 <- as.matrix(exprs(g5)); pd5 <- pData(g5)
storage.mode(e5) <- "numeric"
title5 <- as.character(pd5$title)
parts <- strsplit(title5, "-")
dis5 <- sapply(parts, `[`, 1)        # OA / Normal
reg5 <- sapply(parts, `[`, 2)        # MT / LT
don5 <- sapply(parts, `[`, 3)
stopifnot(all(dis5 %in% c("OA", "Normal")))

map5 <- ann[probe %in% rownames(e5)]
mu <- rowMeans(e5[map5$probe, ], na.rm = TRUE)
map5[, mu := mu]
map5 <- map5[order(-mu)]
map5 <- map5[!duplicated(symbol)]
eg5 <- e5[map5$probe, , drop = FALSE]
rownames(eg5) <- map5$symbol
eg5 <- eg5[is.finite(rowSums(eg5)), ]
cat("GSE51588 基因级矩阵:", dim(eg5), "\n")

# ---- R1: 区域配对 (仅 OA, 供体阻断) ----
oa <- which(dis5 == "OA")
X1 <- model.matrix(~ factor(don5[oa]) + factor(reg5[oa], levels = c("LT", "MT")))
fit1 <- eBayes(lmFit(eg5[, oa], X1))
tt1 <- topTable(fit1, coef = "factor(reg5[oa], levels = c(\"LT\", \"MT\"))MT",
                number = Inf, sort.by = "none")
r1 <- setNames(tt1$t, rownames(tt1))
cat("R1 区域排序: 上调(MT) ", sum(r1 > 0), " / 下调 ", sum(r1 < 0), "\n")
fg1 <- run_fgsea(r1, "R1_GSE51588_Medial_vs_Lateral")

# ---- R2: 疾病主效应 (校正区域) ----
X2 <- model.matrix(~ factor(reg5) + factor(dis5, levels = c("Normal", "OA")))
fit2 <- eBayes(lmFit(eg5, X2))
tt2 <- topTable(fit2, coef = "factor(dis5, levels = c(\"Normal\", \"OA\"))OA",
                number = Inf, sort.by = "none")
r2 <- setNames(tt2$t, rownames(tt2))
fg2 <- run_fgsea(r2, "R2_GSE51588_OA_vs_Normal")

# ---- R5: PIEZO1 共表达排序 (OA 样本 Spearman) ----
pzv <- eg5["PIEZO1", oa]
rr <- apply(eg5[, oa], 1, function(x) suppressWarnings(cor(pzv, x, method = "spearman")))
r5 <- rr[!is.na(rr)]
fg5 <- run_fgsea(r5, "R5_GSE51588_PIEZO1_coexpr")

## ============================================================
## ---------- 滑膜数据集 (Affymetrix) ----------
## ============================================================
load_syno <- function(f, platform) {
  g <- getGEO(filename = f, getGPL = FALSE)
  e <- as.matrix(exprs(g)); pd <- pData(g)
  storage.mode(e) <- "numeric"
  if (max(e, na.rm = TRUE) > 100) { e[e < 1] <- 1; e <- log2(e) }
  sym <- if (platform == "GPL96") toTable(hgu133aSYMBOL) else toTable(hgu133plus2SYMBOL)
  map <- data.table(probe = sym$probe_id, symbol = sym$symbol)
  map <- map[probe %in% rownames(e) & symbol != ""]
  mu <- rowMeans(e[map$probe, ], na.rm = TRUE)
  map[, mu := mu]
  map <- map[order(-mu)][!duplicated(symbol)]
  eg <- e[map$probe, , drop = FALSE]; rownames(eg) <- map$symbol
  # 鲁棒分组: 扫描 characteristics 列, 剔除 RA
  grp <- rep(NA, nrow(pd))
  for (cc in colnames(pd)) {
    v <- as.character(pd[[cc]])
    if (any(grepl("osteoarthrit", v, ignore.case = TRUE))) {
      vv <- toupper(trimws(sub(".*:", "", v)))
      g2 <- ifelse(grepl("OSTEOARTH", vv) | vv == "OA", "OA",
                   ifelse(grepl("HEALTHY", vv) | grepl("NORMAL", vv), "Healthy", NA))
      g2[grepl("RHEUM", vv)] <- NA
      if (sum(!is.na(g2)) >= 10) { grp <- g2; break }
    }
  }
  list(eg = eg, grp = grp)
}

# ---- R3: GSE55235 ----
s3 <- load_syno("data/GSE55235_series_matrix.txt.gz", "GPL96")
keep <- !is.na(s3$grp)
X3 <- model.matrix(~ factor(s3$grp[keep], levels = c("Healthy", "OA")))
fit3 <- eBayes(lmFit(s3$eg[, keep], X3))
tt3 <- topTable(fit3, coef = 2, number = Inf, sort.by = "none")
r3 <- setNames(tt3$t, rownames(tt3))
cat("GSE55235: OA", sum(s3$grp == "OA"), "Healthy", sum(s3$grp == "Healthy"), "\n")
fg3 <- run_fgsea(r3, "R3_GSE55235_OA_vs_Healthy")

# ---- R4: GSE82107 ----
s4 <- load_syno("data/GSE82107_series_matrix.txt.gz", "GPL570")
keep4 <- !is.na(s4$grp)
X4 <- model.matrix(~ factor(s4$grp[keep4], levels = c("Healthy", "OA")))
fit4 <- eBayes(lmFit(s4$eg[, keep4], X4))
tt4 <- topTable(fit4, coef = 2, number = Inf, sort.by = "none")
r4 <- setNames(tt4$t, rownames(tt4))
cat("GSE82107: OA", sum(s4$grp == "OA"), "Healthy", sum(s4$grp == "Healthy"), "\n")
fg4 <- run_fgsea(r4, "R4_GSE82107_OA_vs_Healthy")

## ---------- 汇总 ----------
allfg <- rbindlist(list(fg1, fg2, fg3, fg4, fg5), fill = TRUE)
allfg[, comparison := factor(comparison, levels = c(
  "R1_GSE51588_Medial_vs_Lateral", "R2_GSE51588_OA_vs_Normal",
  "R3_GSE55235_OA_vs_Healthy", "R4_GSE82107_OA_vs_Healthy",
  "R5_GSE51588_PIEZO1_coexpr"))]
fwrite(allfg[, .(comparison, pathway, NES, size, pval, padj, leadingEdge =
                   sapply(leadingEdge, paste, collapse = ";"))],
       file.path(resd, "09_GSEA_Hallmark_all.csv"))

cat("\n===== 各比较 FDR<0.05 Top 通路 =====\n")
for (cc in levels(allfg$comparison)) {
  t <- allfg[comparison == cc & padj < 0.05][order(-abs(NES))]
  cat("\n--", cc, "--\n")
  if (nrow(t)) print(t[1:8, .(pathway, NES, padj)]) else cat("(无 FDR<0.05)\n")
}

## ============================================================
## ---------- 图 ----------
## ============================================================
## fig30: R1 区域配对气泡图
g30 <- bubble(fg1, "Subchondral bone Medial vs Lateral GSEA (Hallmark)",
              subtitle = "GSE51588 | 20 OA donor-paired limma rank | NES>0 = enriched in weight-bearing region")
ggsave(file.path(figd, "fig30_GSEA_region_paired.png"), g30, width = 10, height = 7, dpi = 300)

## fig31: R2 疾病主效应
g31 <- bubble(fg2, "Subchondral bone OA vs Normal GSEA (Hallmark)",
              subtitle = "GSE51588 | adjusted region main effect | NES>0 = OA enrichment")
ggsave(file.path(figd, "fig31_GSEA_disease.png"), g31, width = 10, height = 7, dpi = 300)

## fig32: R3+R4 滑膜并排
fg34 <- rbindlist(list(fg3, fg4))
lab <- c(R3_GSE55235_OA_vs_Healthy = "GSE55235 synovium",
         R4_GSE82107_OA_vs_Healthy = "GSE82107 synovium")
d32 <- as.data.table(fg34)[padj < 0.25][order(padj)]
d32[, laby := gsub("HALLMARK_|_", " ", gsub("HALLMARK_", "", pathway))]
d32[, laby := gsub("_", " ", laby)]
d32$ds <- lab[as.character(d32$comparison)]
topp <- unique(d32[order(padj)]$laby)[1:14]
d32 <- d32[laby %in% topp]
d32[, laby := factor(laby, levels = sort(unique(laby)))]
g32 <- ggplot(d32, aes(NES, laby, color = padj)) +
  geom_point(aes(size = size), alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  scale_color_gradientn(colors = c("#B2182B", "#EF8A62", "#FDDBC7", "#D1E5F0", "#67A9CF", "#2166AC"),
                         trans = "reverse", name = "FDR q") +
  scale_size_continuous(range = c(2, 7), name = "Gene count") +
  facet_grid(ds ~ ., scales = "free_y", space = "free_y") +
  labs(x = "NES (right = OA enrichment)", y = NULL, title = "Synovium OA vs Healthy GSEA (Hallmark)") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 8.5))
ggsave(file.path(figd, "fig32_GSEA_synovium.png"), g32, width = 10, height = 8, dpi = 300)

## fig33: R5 PIEZO1 共表达
g33 <- bubble(fg5, "PIEZO1 co-expression program GSEA (Hallmark)",
              subtitle = "GSE51588 | 40 OA subchondral bone samples | Spearman rank | NES>0 = increases with PIEZO1")
ggsave(file.path(figd, "fig33_GSEA_piezo1_coexpr.png"), g33, width = 10, height = 7, dpi = 300)

## fig34: NES 跨比较汇总热图 (pheatmap)
suppressMessages(library(pheatmap))
sigp <- allfg[padj < 0.25, .(n = .N), by = pathway][n >= 2][order(-n)]$pathway
if (length(sigp) > 25) sigp <- allfg[pathway %in% sigp][, .(m = mean(padj)), by = pathway][order(m)][1:25]$pathway
nes <- dcast(allfg[pathway %in% sigp], pathway ~ comparison, value.var = "NES")
mat <- as.matrix(nes[, -1]); rownames(mat) <- gsub("HALLMARK_", "", nes$pathway)
colnames(mat) <- c("Subchondral\nMedial vs Lateral", "Subchondral\nOA vs Normal",
                   "Synovium55235\nOA vs Heal", "Synovium82107\nOA vs Heal", "PIEZO1\nCo-expression")
bk <- seq(-3, 3, length.out = 61)
png(file.path(figd, "fig34_GSEA_NES_heatmap.png"), width = 2200, height = 2400, res = 300)
pheatmap(mat, breaks = bk,
         color = colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(60),
         border_color = NA, cluster_cols = FALSE,
         main = "Hallmark pathway NES across-comparison heatmap", fontsize_row = 8.5, fontsize = 10,
         angle_col = 45)
dev.off()

## fig35: 关键通路 GSEA 运行分数曲线 (PIEZO1 共表达 + 区域配对)
keyp <- intersect(c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                    "HALLMARK_INFLAMMATORY_RESPONSE",
                    "HALLMARK_IL6_JAK_STAT3_SIGNALING"), names(paths))
plots <- list()
for (kp in keyp) {
  p1 <- plotEnrichment(stats = r5, pathway = paths[[kp]]) +
    labs(title = NULL, x = NULL) + theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 9))
  plots[[kp]] <- ggplotGrob(p1)
}
library(gridExtra)
g35 <- arrangeGrob(grobs = plots, ncol = 2,
                   top = "PIEZO1 co-expression ranking - GSEA running score of key inflammatory pathways (GSE51588 OA)")
ggsave(file.path(figd, "fig35_GSEA_running_score.png"), g35, width = 10, height = 6, dpi = 300)

cat("\n全部 GSEA 图表输出完毕\n")
