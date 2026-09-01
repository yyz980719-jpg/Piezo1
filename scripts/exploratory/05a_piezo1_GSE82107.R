# ============================================================
#  05a_piezo1_GSE82107.R — PIEZO1 快查: 第二个滑膜数据集
#  GSE82107: 滑膜, 7 Healthy vs 10 OA, GPL570 (HG-U133 Plus 2.0)
#  输出: figs/GSE82107_PIEZO1_boxplot.png / _ROC.png
#  复用方法: 改文件名与分组标签即可用于其他 GPL570 数据集
# ============================================================
suppressMessages({
  library(GEOquery)
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(ggplot2)
  library(ggpubr)
  library(pROC)
})
source("scripts/_bootstrap.R")
d <- PIEZO1_ROOT
figs <- file.path(d, "figs")
targets <- c("PIEZO1", "FAM38A")

## ---------- 读入 ----------
g <- getGEO(filename = file.path(d, "data", "GSE82107_series_matrix.txt.gz"), getGPL = FALSE)
e <- as.matrix(exprs(g)); pd <- pData(g)
cat(sprintf("表达矩阵: %d 探针 x %d 样本\n", nrow(e), ncol(e)))

## ---------- 分组（鲁棒提取：扫描含 OA 值的列）----------
grp_raw <- NULL
for (cc in colnames(pd)) {
  v <- as.character(pd[[cc]])
  if (any(grepl("osteoarthrit", v, ignore.case = TRUE))) { grp_raw <- v; break }
}
vv <- toupper(trimws(sub(".*:", "", grp_raw)))
grp <- ifelse(vv == "OA" | grepl("OSTEOARTH", vv), "OA", "Healthy")
names(grp) <- rownames(pd)
cat(sprintf("分组: %s\n", paste(names(table(grp)), table(grp), sep = "=", collapse = ", ")))
if (max(e, na.rm = TRUE) > 100) { e <- log2(e + 1); cat("已做 log2 变换\n") }

## ---------- PIEZO1 探针 ----------
pp <- intersect(subset(suppressWarnings(toTable(hgu133plus2SYMBOL)),
                        symbol %in% targets)$probe_id, rownames(e))
cat("PIEZO1 探针:", paste(pp, collapse = ", "), "\n")
for (pb in pp) {
  v1 <- as.numeric(e[pb, ])
  t2 <- t.test(v1[grp == "OA"], v1[grp == "Healthy"])
  cat(sprintf("  探针 %s: OA %.2f vs Healthy %.2f | 差值 %+.2f | P = %.3g\n",
              pb, mean(v1[grp == "OA"]), mean(v1[grp == "Healthy"]),
              mean(v1[grp == "OA"]) - mean(v1[grp == "Healthy"]), t2$p.value))
}

## ---------- 作图（取平均表达最高的探针）----------
best <- pp[which.max(rowMeans(e[pp, , drop = FALSE]))]
cat("作图探针:", best, "\n")
df <- data.frame(group = factor(grp, levels = c("Healthy", "OA")),
                 expr = as.numeric(e[best, ]))
tt <- t.test(expr ~ group, df)
m <- tapply(df$expr, df$group, mean)
cat(sprintf("主探针 %s: Healthy %.2f vs OA %.2f | P = %.3g\n",
            best, m["Healthy"], m["OA"], tt$p.value))

p <- ggboxplot(df, x = "group", y = "expr", color = "group", add = "jitter",
               palette = c("#2f6eb5", "#d64541"), xlab = NULL,
               ylab = "PIEZO1 expression",
               title = paste0("PIEZO1 in OA synovium - GSE82107 (", best, ")")) +
  stat_compare_means(method = "t.test", size = 3.2)
ggsave(file.path(figs, "GSE82107_PIEZO1_boxplot.png"), p, width = 4.6, height = 4.6, dpi = 300)
cat("图已保存: GSE82107_PIEZO1_boxplot.png\n")

roc1 <- roc(response = as.integer(df$group == "OA"), predictor = df$expr,
            quiet = TRUE, direction = "auto")
a <- as.numeric(auc(roc1))
cat(sprintf("ROC: AUC = %.3f\n", a))
png(file.path(figs, "GSE82107_PIEZO1_ROC.png"), width = 1600, height = 1600, res = 300)
par(mar = c(4.2, 4.2, 3, 1))
plot(roc1, col = "#d64541", lwd = 2.5, legacy.axes = TRUE,
     main = sprintf("PIEZO1 ROC - GSE82107 synovium\nAUC = %.3f", a),
     xlab = "1 - Specificity", ylab = "Sensitivity")
abline(0, 1, lty = 2, col = "grey60")
dev.off()
cat("图已保存: GSE82107_PIEZO1_ROC.png\n")
