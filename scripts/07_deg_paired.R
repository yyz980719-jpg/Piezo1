## ============================================================
## 07_deg_paired.R —— GSE51588 软骨下骨全转录组差异分析（论文 Figure 2 核心）
##   对比1[主] : OA 内 负重(Medial) vs 非负重(Lateral)，供体配对
##   对比2     : OA vs Normal（校正区域）—— 疾病主效应
##   对比3     : Δ(负重−非负重) 的 OA vs Normal —— 机械响应交互效应
##   模块      : 与 PIEZO1 Δ 共变的基因（Spearman）+ GO/KEGG 富集
## ============================================================
suppressMessages({
  library(GEOquery); library(limma); library(data.table)
  library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db)
})
source("scripts/_bootstrap.R")

## ---------- 1. 数据与样本信息 ----------
g <- getGEO(filename = "data/GSE51588_series_matrix.txt.gz", getGPL = FALSE)
e <- as.matrix(exprs(g)); pd <- pData(g)
ttl <- as.character(pd$title)
dis <- ifelse(grepl("^OA-", ttl), "OA", "Normal")
reg <- ifelse(grepl("-MT-", ttl), "Medial", "Lateral")
don <- sub("^[^-]+-[^-]+-", "", ttl)
uid <- paste(dis, don, sep = "_")
cat(sprintf("[样本] %d | OA %d / Normal %d | Medial %d / Lateral %d | 供体 %d\n",
            ncol(e), sum(dis == "OA"), sum(dis == "Normal"),
            sum(reg == "Medial"), sum(reg == "Lateral"), length(unique(uid))))
if (max(e, na.rm = TRUE) > 100) { e <- log2(e + 1); cat("[尺度] 已 log2 转换\n")
} else cat("[尺度] 数据已为 log2 尺度\n")

## ---------- 2. 注释 → 基因级矩阵（每基因取均值最高探针） ----------
ann <- fread(data_path("GPL13497_probe2symbol.tsv"), header = TRUE,
             sep = "\t", colClasses = "character")
ann <- ann[!is.na(symbol) & nzchar(symbol) &
           grepl("^[A-Za-z][A-Za-z0-9@\\.\\-]*$", symbol)]
common <- intersect(rownames(e), ann$probe)
eg <- e[common, , drop = FALSE]
gl <- ann$symbol[match(common, ann$probe)]
me <- rowMeans(eg, na.rm = TRUE)
o <- order(-me); eg <- eg[o, , drop = FALSE]; gl <- gl[o]
kd <- !duplicated(gl)
emat <- eg[kd, , drop = FALSE]; rownames(emat) <- gl[kd]
storage.mode(emat) <- "numeric"
cat(sprintf("[注释] 探针 %d → 基因 %d | 含 PIEZO1: %s\n",
            length(common), nrow(emat), "PIEZO1" %in% rownames(emat)))

## ---------- 3. 对比1（主分析）：OA 内配对 负重 vs 非负重 ----------
oa <- which(dis == "OA")
dsg <- model.matrix(~ factor(uid[oa]) + factor(reg[oa], levels = c("Lateral", "Medial")))
colnames(dsg)[ncol(dsg)] <- "Medial"
fit <- eBayes(lmFit(emat[, oa], dsg), trend = TRUE)
tt <- topTable(fit, coef = "Medial", number = Inf, sort.by = "P")
tt$gene <- rownames(tt)
fwrite(tt, "results/07_DEG_paired_MedialVsLateral_OA.csv")
n1 <- sum(tt$adj.P.Val < 0.05)
n1u <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0)
cat(sprintf("\n[对比1 OA配对 负重vs非负重] FDR<0.05: %d 个基因 (负重区升高 %d / 降低 %d)\n",
            n1, n1u, n1 - n1u))
i <- which(rownames(tt) == "PIEZO1")
cat(sprintf("  PIEZO1: logFC=%+.3f  P=%.3g  FDR=%.3g  按P排名 %d/%d\n",
            tt$logFC[i], tt$P.Value[i], tt$adj.P.Val[i], i, nrow(tt)))
cat("  Top10 基因:\n")
print(head(tt[, c("gene", "logFC", "P.Value", "adj.P.Val")], 10), row.names = FALSE)

## ---------- 4. 通用火山图 ----------
volcano_plot <- function(tb, title, uplab, dnlab, lfc_cut = 0.25, fdr_cut = 0.05, nlab = 12) {
  d <- data.frame(logFC = tb$logFC, P = tb$P.Value, FDR = tb$adj.P.Val, gene = rownames(tb))
  d$g <- ifelse(d$FDR < fdr_cut & d$logFC > lfc_cut, uplab,
         ifelse(d$FDR < fdr_cut & d$logFC < -lfc_cut, dnlab, "ns"))
  cols <- c("#B2182B", "#2166AC", "grey70"); names(cols) <- c(uplab, dnlab, "ns")
  sig <- d[d$g != "ns", ]
  key <- unique(rbind(head(sig[order(-abs(sig$logFC)), ], nlab), d[d$gene == "PIEZO1", ]))
  p <- ggplot(d, aes(logFC, -log10(P))) +
    geom_point(aes(color = g), size = 0.6, alpha = 0.5) +
    scale_color_manual(values = cols, drop = FALSE) +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 2, color = "grey60") +
    labs(title = title, x = "log2FC", y = "-log10(P)", color = NULL) +
    theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 13))
  if (nrow(key) > 0) {
    p <- p + ggrepel::geom_text_repel(data = key, aes(label = gene),
             size = 3, max.overlaps = 25, segment.size = 0.3, color = "grey15")
  }
  p
}
fig11 <- volcano_plot(tt, "OA subchondral bone: Medial vs Lateral (donor-paired, n=20 pairs)",
                      "Up in Medial", "Up in Lateral")
ggsave("figs/fig11_volcano_paired.png", fig11, width = 8, height = 6.5, dpi = 160)

## ---------- 5. Top50 热图（全部样本） ----------
topg <- tt$gene[tt$adj.P.Val < 0.05]
if (length(topg) < 30) topg <- head(tt$gene, 50) else topg <- head(topg, 50)
M <- emat[topg, , drop = FALSE]
sord <- order(dis, reg, uid)
M <- M[, sord]
andf <- data.frame(Disease = dis[sord], Region = reg[sord], row.names = colnames(M))
pheatmap::pheatmap(t(scale(t(M))), annotation_col = andf, show_colnames = FALSE,
  fontsize_row = 7, border_color = NA,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = sprintf("Paired DEG analysis Top%d genes (z-score)", length(topg)),
  filename = "figs/fig12_heatmap_top50.png", width = 10, height = 9)

## ---------- 6. Top6 DEG 配对折线图 ----------
idx_m <- which(dis == "OA" & reg == "Medial"); idx_l <- which(dis == "OA" & reg == "Lateral")
pairs_oa <- intersect(uid[idx_m], uid[idx_l])
mp <- idx_m[match(pairs_oa, uid[idx_m])]; lp <- idx_l[match(pairs_oa, uid[idx_l])]
top6 <- head(tt$gene, 6)
png("figs/fig13_topDEG_pairedlines.png", width = 1800, height = 1500, res = 160)
par(mfrow = c(3, 2), mar = c(4, 4.2, 2.5, 1))
for (gn in top6) {
  mv <- emat[gn, mp]; lv <- emat[gn, lp]
  pr <- tt$P.Value[tt$gene == gn][1]
  plot(1:2, range(c(mv, lv)), type = "n", xaxt = "n", xlab = "",
       ylab = "log2 expression", xlim = c(0.8, 2.2),
       main = sprintf("%s   P=%.2g", gn, pr), cex.main = 1.0)
  for (k in seq_along(mv)) segments(1, lv[k], 2, mv[k], col = grey(0.75), lwd = 1)
  points(rep(1, length(lv)), lv, pch = 16, cex = 1.1)
  points(rep(2, length(mv)), mv, pch = 16, cex = 1.1, col = "#B2182B")
  axis(1, at = 1:2, labels = c("Non-weight-bearing\nLateral", "Weight-bearing\nMedial"))
}
dev.off()

## ---------- 7. GO / KEGG 富集 ----------
pick_genes <- function(tb, direction, min_n = 30) {
  s <- sign(tb$logFC) == direction
  for (cut in c(0.25, 0.15, 0.05)) {
    g <- tb$gene[tb$adj.P.Val < 0.05 & s & tb$logFC * direction > cut]
    if (length(g) >= min_n) return(list(genes = g, cut = cut))
  }
  list(genes = head(tb$gene[tb$adj.P.Val < 0.05 & s], min_n), cut = 0)
}
run_go <- function(genes, file, title, figfile) {
  if (length(genes) < 10) { cat(sprintf("  [GO] %s: 基因数 %d <10，跳过\n", title, length(genes))); return(invisible(NULL)) }
  ego <- suppressWarnings(enrichGO(gene = genes, keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
      ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1))
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    fwrite(as.data.frame(ego), file)
    p <- dotplot(ego, showCategory = 15) + ggtitle(title) +
      theme(axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold", size = 13))
    ggsave(figfile, p, width = 9, height = 7.5, dpi = 160)
    cat(sprintf("  [GO] %s: %d 条通路 → %s\n", title, nrow(as.data.frame(ego)), figfile))
    print(head(as.data.frame(ego)[, c("Description", "pvalue", "p.adjust", "Count")], 5))
  } else cat(sprintf("  [GO] %s: 无显著通路\n", title))
  invisible(ego)
}
run_kegg <- function(genes, file, figfile) {
  if (length(genes) < 10) return(invisible(NULL))
  egids <- tryCatch(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
      OrgDb = org.Hs.eg.db)$ENTREZID, error = function(e) NULL)
  if (length(egids) < 10) { cat("  [KEGG] 映射基因不足，跳过\n"); return(invisible(NULL)) }
  kk <- tryCatch(enrichKEGG(gene = egids, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.1),
                 error = function(e) { cat("  [KEGG] 失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(kk) || nrow(as.data.frame(kk)) == 0) { cat("  [KEGG] 无显著通路\n"); return(invisible(NULL)) }
  fwrite(as.data.frame(kk), file)
  p <- dotplot(kk, showCategory = 15) + ggtitle("DEG KEGG pathways") +
    theme(axis.text.y = element_text(size = 8), plot.title = element_text(face = "bold", size = 13))
  ggsave(figfile, p, width = 9, height = 7, dpi = 160)
  cat(sprintf("  [KEGG] %d 条通路 → %s\n", nrow(as.data.frame(kk)), figfile))
  print(head(as.data.frame(kk)[, c("Description", "pvalue", "p.adjust", "Count")], 8))
  invisible(kk)
}
pu <- pick_genes(tt, +1); dn1 <- pick_genes(tt, -1)
cat(sprintf("\n[富集输入] 负重区升高 %d 基因(logFC>%.2f) | 非负重区升高 %d 基因(logFC<-%.2f)\n",
            length(pu$genes), pu$cut, length(dn1$genes), dn1$cut))
run_go(pu$genes, "results/07_GO_MedialUp.csv", "OA Medial-up genes · GO Biological Process", "figs/fig14_GO_MedialUp.png")
run_go(dn1$genes, "results/07_GO_LateralUp.csv", "OA Lateral-up genes · GO Biological Process", "figs/fig15_GO_LateralUp.png")
alldeg <- tt$gene[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0.25]
if (length(alldeg) < 10) alldeg <- tt$gene[tt$adj.P.Val < 0.05]
run_kegg(alldeg, "results/07_KEGG_paired.csv", "figs/fig16_KEGG_DEG.png")

## ---------- 8. 对比2: OA vs Normal（校正区域） ----------
dsg2 <- model.matrix(~ factor(dis) + factor(reg))
colnames(dsg2) <- c("(Intercept)", "disOA", "regMedial")
fit2 <- eBayes(lmFit(emat, dsg2), trend = TRUE)
tt2 <- topTable(fit2, coef = "disOA", number = Inf, sort.by = "P")
tt2$gene <- rownames(tt2)
fwrite(tt2, "results/07_DEG_disease_OAvsNormal.csv")
n2 <- sum(tt2$adj.P.Val < 0.05)
cat(sprintf("\n[对比2 OA vs Normal（校正区域）] FDR<0.05: %d 个基因\n", n2))
j <- which(rownames(tt2) == "PIEZO1")
cat(sprintf("  PIEZO1: logFC=%+.3f  P=%.3g  FDR=%.3g\n",
            tt2$logFC[j], tt2$P.Value[j], tt2$adj.P.Val[j]))
cat("  Top10 基因:\n")
print(head(tt2[, c("gene", "logFC", "P.Value", "adj.P.Val")], 10), row.names = FALSE)
fig17 <- volcano_plot(tt2, "OA vs Normal subchondral bone (adjusted region, n=40 vs 10)", "Upregulated", "Downregulated")
ggsave("figs/fig17_volcano_disease.png", fig17, width = 8, height = 6.5, dpi = 160)
disdeg <- tt2$gene[tt2$adj.P.Val < 0.05]
run_go(head(disdeg, 1000), "results/07_GO_disease.csv", "OA vs Normal DEGs · GO Biological Process", "figs/fig17b_GO_disease.png")

## ---------- 9. 对比3: 交互效应 Δ(负重−非负重) 的疾病间差异 ----------
uids <- unique(uid)
dl <- setNames(lapply(uids, function(u) {
  m <- which(uid == u & reg == "Medial"); l <- which(uid == u & reg == "Lateral")
  if (length(m) == 1 && length(l) == 1) emat[, m] - emat[, l] else NULL
}), uids)
dl <- dl[!sapply(dl, is.null)]
D <- do.call(cbind, dl)
dd <- sub("_.*", "", colnames(D))
dsg3 <- model.matrix(~ factor(dd))
colnames(dsg3) <- c("(Intercept)", "disOA")
fit3 <- eBayes(lmFit(D, dsg3), trend = TRUE)
tt3 <- topTable(fit3, coef = "disOA", number = Inf, sort.by = "P")
tt3$gene <- rownames(tt3)
fwrite(tt3, "results/07_DEG_interaction.csv")
n3 <- sum(tt3$adj.P.Val < 0.05)
cat(sprintf("\n[对比3 交互: Δ(负重−非负重) OA vs Normal] FDR<0.05: %d 个基因\n", n3))
k <- which(rownames(tt3) == "PIEZO1")
cat(sprintf("  PIEZO1: Δ差异 logFC=%+.3f  P=%.3g  FDR=%.3g\n",
            tt3$logFC[k], tt3$P.Value[k], tt3$adj.P.Val[k]))
cat("  Top10 基因:\n")
print(head(tt3[, c("gene", "logFC", "P.Value", "adj.P.Val")], 10), row.names = FALSE)
fig18 <- volcano_plot(tt3, "Mechanical-response interaction: Δ(Medial−Lateral) difference, OA vs Normal",
                      "Stronger in OA", "Weaker in OA", lfc_cut = 0.2)
ggsave("figs/fig18_volcano_interaction.png", fig18, width = 8, height = 6.5, dpi = 160)

## ---------- 10. PIEZO1 共变模块（Δ 相关） ----------
pz_d <- D["PIEZO1", ]
rho <- apply(D, 1, function(x) suppressWarnings(cor(x, pz_d, method = "spearman")))
cor_tb <- data.frame(gene = rownames(D), rho = rho)
cor_tb <- cor_tb[order(-cor_tb$rho), ]
fwrite(cor_tb, "results/07_PIEZO1_delta_correlation.csv")
posg <- cor_tb$gene[cor_tb$rho >= 0.6]
if (length(posg) < 30) posg <- head(cor_tb$gene, 300)
cat(sprintf("\n[模块] 与 PIEZO1 Δ Spearman ρ≥0.6: %d 基因；取 %d 个做富集\n",
            sum(cor_tb$rho >= 0.6), length(posg)))
ovp <- intersect(tt$gene[tt$adj.P.Val < 0.05], head(cor_tb$gene, 500))
cat(sprintf("  配对 DEG ∩ 模块 Top500: %d 个重叠\n", length(ovp)))

png("figs/fig19_piezo1_module.png", width = 2000, height = 900, res = 160)
par(mfrow = c(1, 2), mar = c(5, 9, 3, 1))
t15 <- head(cor_tb, 15)
barplot(rev(t15$rho), names.arg = rev(t15$gene), horiz = TRUE, las = 1, cex.names = 0.72,
        col = "#4393C3", xlab = "Spearman rho (vs PIEZO1 delta)",
        main = "Top15 genes co-varying with PIEZO1 mechanical response", border = NA)
par(mar = c(5, 4.5, 3, 1))
g1 <- t15$gene[1]
plot(pz_d, D[g1, ], pch = 16, cex = 1.2,
     col = ifelse(dd == "OA", "#B2182B", "#4393C3"),
     xlab = "PIEZO1  delta(Medial-Lateral)", ylab = sprintf("%s  delta", g1),
     main = sprintf("%s co-varying with PIEZO1 (rho=%.2f)", g1, t15$rho[1]))
abline(lm(D[g1, ] ~ pz_d), lty = 2, col = "grey30")
legend("topleft", c("OA donor", "Normal donor"), col = c("#B2182B", "#4393C3"), pch = 16, bty = "n")
dev.off()
run_go(posg, "results/07_GO_piezo1_module.csv", "PIEZO1 co-varying mechanical-response module · GO Biological Process", "figs/fig20_GO_module.png")

## ---------- 11. PIEZO1 在全基因组配对分析中的排名 ----------
i <- which(rownames(tt) == "PIEZO1")
dpi <- data.frame(logFC = tt$logFC[i], P = tt$P.Value[i])
p <- ggplot(data.frame(logFC = tt$logFC, P = tt$P.Value), aes(logFC, -log10(P))) +
  geom_point(color = "grey65", size = 0.5, alpha = 0.5) +
  geom_point(data = dpi, aes(logFC, -log10(P)), color = "#B2182B", size = 3.5) +
  geom_hline(yintercept = -log10(dpi$P), linetype = 3, color = "#B2182B", alpha = 0.5) +
  geom_vline(xintercept = dpi$logFC, linetype = 3, color = "#B2182B", alpha = 0.5) +
  labs(title = "Position of PIEZO1 in the genome-wide paired analysis",
       subtitle = sprintf("PIEZO1: logFC=%+.2f, P=%.2g, FDR=%.2g | rank by P-value %d / %d (top %.1f%%)",
                          tt$logFC[i], tt$P.Value[i], tt$adj.P.Val[i], i, nrow(tt), 100 * i / nrow(tt)),
       x = "log2FC (Medial - Lateral)", y = "-log10(P)") +
  theme_bw(base_size = 12) +
  ggrepel::geom_text_repel(data = data.frame(logFC = dpi$logFC, P = dpi$P, gene = "PIEZO1"),
      aes(logFC, -log10(P), label = gene), color = "#B2182B", size = 4, fontface = 2)
ggsave("figs/fig21_piezo1_rank.png", p, width = 8, height = 6, dpi = 160)

## ---------- 12. 汇总 ----------
cat("\n================= 汇 总 =================\n")
cat(sprintf("对比1 OA配对 负重vs非负重 : %d DEG（FDR<0.05，其中负重区升高 %d）\n", n1, n1u))
cat(sprintf("对比2 OA vs Normal         : %d DEG\n", n2))
cat(sprintf("对比3 机械响应交互         : %d DEG\n", n3))
cat(sprintf("PIEZO1 共变模块（ρ≥0.6）    : %d 基因；与配对 DEG 重叠 %d 个\n",
            sum(cor_tb$rho >= 0.6), length(ovp)))
cat("输出: figs/fig11~fig21 + results/07_*.csv\n")
cat("=========================================\n")
