# ============================================================
#  05b_piezo1_GSE51588.R — PIEZO1 快查: 软骨下骨数据集
#  GSE51588: 胫骨平台软骨下骨, 10 Normal vs 40 OA, Agilent GPL13497
#  注释: GPL13497_family.soft.gz 平台表 (ID -> 基因符号)
#  设计说明: 20 例 OA 供体 + 5 例非 OA 供体, 各取内外侧平台
#           (OA n=40, Normal n=10; 同一供体贡献 2 个区域样本)
# ============================================================
suppressMessages({
  library(GEOquery)
  library(ggplot2)
  library(ggpubr)
  library(pROC)
  library(data.table)
})
source("scripts/_bootstrap.R")
d <- PIEZO1_ROOT
figs <- file.path(d, "figs")
targets <- c("PIEZO1", "FAM38A")

## ---------- 1. 读取 Agilent 平台注释 ----------
## 完整 SOFT 文件下载极慢且曾截断; PIEZO1 行已验证完整,
## 用迷你注释表 (探针->符号, 从 SOFT 提取: 表头 + PIEZO1/FAM38A 行)
mini <- file.path(d, "data", "GPL13497_probe2symbol.tsv")
ann <- fread(mini, header = TRUE)
cat("迷你注释列名:", paste(colnames(ann), collapse = ", "), "\n")
pz <- unique(ann[[1]][ann[[2]] %in% targets])
cat("PIEZO1 探针:", paste(pz, collapse = ", "), "\n")

## ---------- 2. 读入表达数据 ----------
g <- getGEO(filename = file.path(d, "data", "GSE51588_series_matrix.txt.gz"), getGPL = FALSE)
e <- as.matrix(exprs(g)); pd <- pData(g)
cat(sprintf("表达矩阵: %d 探针 x %d 样本\n", nrow(e), ncol(e)))

## ---------- 3. 分组 ----------
grp_raw <- NULL
for (cc in colnames(pd)) {
  v <- as.character(pd[[cc]])
  if (any(grepl("disease state", v))) { grp_raw <- v; break }
}
vv <- toupper(trimws(sub(".*:", "", grp_raw)))
grp <- ifelse(vv == "OA" | grepl("OSTEOARTH", vv), "OA", "Normal")
names(grp) <- rownames(pd)
cat(sprintf("分组: %s\n", paste(names(table(grp)), table(grp), sep = "=", collapse = ", ")))
if (max(e, na.rm = TRUE) > 100) { e <- log2(e + 1); cat("已做 log2 变换\n") }

## ---------- 4. PIEZO1 检验 ----------
pp <- intersect(pz, rownames(e))
stopifnot(length(pp) > 0)
for (pb in pp) {
  v1 <- as.numeric(e[pb, ])
  t2 <- t.test(v1[grp == "OA"], v1[grp == "Normal"])
  cat(sprintf("  探针 %s: OA %.2f vs Normal %.2f | 差值 %+.2f | P = %.3g\n",
              pb, mean(v1[grp == "OA"]), mean(v1[grp == "Normal"]),
              mean(v1[grp == "OA"]) - mean(v1[grp == "Normal"]), t2$p.value))
}

## ---------- 5. 作图（平均表达最高的探针）----------
best <- pp[which.max(rowMeans(e[pp, , drop = FALSE]))]
cat("作图探针:", best, "\n")
df <- data.frame(group = factor(grp, levels = c("Normal", "OA")),
                 expr = as.numeric(e[best, ]))
tt <- t.test(expr ~ group, df)
m <- tapply(df$expr, df$group, mean)
cat(sprintf("主探针 %s: Normal %.2f vs OA %.2f | P = %.3g\n",
            best, m["Normal"], m["OA"], tt$p.value))

p <- ggboxplot(df, x = "group", y = "expr", color = "group", add = "jitter",
               palette = c("#2f6eb5", "#d64541"), xlab = NULL,
               ylab = "PIEZO1 expression (normalized)",
               title = paste0("PIEZO1 in OA subchondral bone - GSE51588 (", best, ")")) +
  stat_compare_means(method = "t.test", size = 3.2)
ggsave(file.path(figs, "GSE51588_PIEZO1_boxplot.png"), p, width = 4.6, height = 4.6, dpi = 300)
cat("图已保存: GSE51588_PIEZO1_boxplot.png\n")

roc1 <- roc(response = as.integer(df$group == "OA"), predictor = df$expr,
            quiet = TRUE, direction = "auto")
a <- as.numeric(auc(roc1))
cat(sprintf("ROC: AUC = %.3f\n", a))
png(file.path(figs, "GSE51588_PIEZO1_ROC.png"), width = 1600, height = 1600, res = 300)
par(mar = c(4.2, 4.2, 3, 1))
plot(roc1, col = "#d64541", lwd = 2.5, legacy.axes = TRUE,
     main = sprintf("PIEZO1 ROC - GSE51588 subchondral bone\nAUC = %.3f", a),
     xlab = "1 - Specificity", ylab = "Sensitivity")
abline(0, 1, lty = 2, col = "grey60")
dev.off()
cat("图已保存: GSE51588_PIEZO1_ROC.png\n")
