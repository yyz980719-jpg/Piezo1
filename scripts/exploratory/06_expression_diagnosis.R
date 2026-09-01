# ============================================================
# 06_expression_diagnosis.R — PIEZO1 表达与诊断层 完整分析
# 数据:
#   GSE51588 软骨下骨 (5 Normal 供体 + 20 OA 供体; 每供体
#            内侧平台 MT=负重区 / 外侧平台 LT=非负重区, 配对设计)
#   GSE82107 滑膜   (7 Healthy + 10 OA, GPL570, 3 个 PIEZO1 探针)
#   GSE55235 滑膜   (10 Healthy + 10 OA, GPL96, 探针 202771_at)
# 分析角度:
#   A1 疾病差异      OA vs Normal/Healthy (三数据集)
#   A2 机械应力      同一供体 负重Medial vs 非负重Lateral 配对 (GSE51588 核心)
#   A3 疾病x应力交互 delta(MT-LT) 在 OA vs Normal 间比较
#   A4 年龄相关      Spearman (GSE51588)
#   A5 诊断效能      ROC / AUC (三数据集)
#   A6 多探针一致性  GSE82107 三探针
#   A7 跨数据集森林图 + 汇总表
# 输出: figs/fig01..fig10 + results/PIEZO1_expression_diagnosis_summary.csv
# ============================================================
suppressMessages({
  library(GEOquery); library(ggplot2); library(ggpubr)
  library(pROC); library(data.table)
})
source("scripts/_bootstrap.R")
d    <- PIEZO1_ROOT
figs <- file.path(d, "figs"); res <- file.path(d, "results")
dir.create(figs, recursive = TRUE, showWarnings = FALSE)
dir.create(res,  recursive = TRUE, showWarnings = FALSE)
BLUE <- "#2f6eb5"; RED <- "#d64541"; TEAL <- "#0e8f6e"; PURP <- "#7b52a1"
sumrows <- list()   # 汇总表收集器
addrow <- function(...) sumrows[[length(sumrows)+1]] <<- data.frame(...)

# ============================================================
# PART A: GSE51588 软骨下骨 (核心: 疾病 + 机械应力配对)
# ============================================================
cat("========== A. GSE51588 软骨下骨 ==========\n")
g <- getGEO(filename = file.path(d, "data", "GSE51588_series_matrix.txt.gz"),
            getGPL = FALSE)
e <- as.matrix(exprs(g)); pd <- pData(g)
if (max(e, na.rm = TRUE) > 100) e <- log2(e + 1)
title <- as.character(pd$title)
dis   <- ifelse(grepl("^OA-", title), "OA", "Normal")
reg   <- ifelse(grepl("-MT-", title), "Medial", "Lateral")
donum <- sub("^(Normal|OA)-(LT|MT)-", "", title)
donor <- paste(dis, donum, sep = "_")
stopifnot(all(table(donor) == 2))          # 每个供体恰好 2 个部位样本
age   <- suppressWarnings(as.numeric(as.character(pd[, "age:ch1"])))
pz    <- "A_23_P140738"                    # PIEZO1 唯一探针 (GPL13497)
v     <- as.numeric(e[pz, ]); names(v) <- rownames(pd)
g4    <- factor(paste(dis, reg, sep = "\n"))   # 四分组
cat(sprintf("样本: OA %d (供体 %d), Normal %d (供体 %d); 探针 %s\n",
            sum(dis == "OA"), length(unique(donor[dis == "OA"])),
            sum(dis == "Normal"), length(unique(donor[dis == "Normal"])), pz))

## ---- A1 疾病差异: OA vs Normal ----
t1 <- t.test(v[dis == "OA"], v[dis == "Normal"])
d1 <- mean(v[dis == "OA"]) - mean(v[dis == "Normal"])   # OA - Normal (方向与CI一致)
roc1 <- roc(response = as.integer(dis == "OA"), predictor = v,
            quiet = TRUE, direction = "auto"); auc1 <- as.numeric(auc(roc1))
cat(sprintf("[A1] OA vs Normal: %+.3f (95%%CI %+.3f~%+.3f) P=%.3g | AUC=%.3f\n",
            d1, t1$conf.int[1], t1$conf.int[2], t1$p.value, auc1))
addrow(Dataset = "GSE51588", Tissue = "Subchondral bone", Comparison = "OA vs Normal",
       Probe = pz, n_Ctrl = sum(dis == "Normal"), n_Case = sum(dis == "OA"),
       Diff = d1, CI_low = t1$conf.int[1], CI_high = t1$conf.int[2],
       P = t1$p.value, AUC = auc1, Method = "Welch t-test",
       Direction = "OA - Normal")

## ---- A2 机械应力: 同一供体 负重Medial vs 非负重Lateral (配对) ----
w <- data.frame(donor = donor, dis = dis, reg = reg, expr = v)
m <- merge(w[w$reg == "Medial", c("donor", "dis", "expr")],
           w[w$reg == "Lateral", c("donor", "expr")], by = "donor",
           suffixes = c("_M", "_L"))
m$delta <- m$expr_M - m$expr_L
tp_all <- t.test(m$expr_M, m$expr_L, paired = TRUE)
tp_oa  <- t.test(m$expr_M[m$dis == "OA"], m$expr_L[m$dis == "OA"], paired = TRUE)
tp_no  <- t.test(m$expr_M[m$dis == "Normal"], m$expr_L[m$dis == "Normal"], paired = TRUE)
cat(sprintf("[A2] Medial vs Lateral 全部配对(n=%d): %+.3f (95%%CI %+.3f~%+.3f) P=%.3g\n",
            nrow(m), tp_all$estimate, tp_all$conf.int[1], tp_all$conf.int[2], tp_all$p.value))
cat(sprintf("[A2] Medial vs Lateral 仅OA患者(n=%d对): %+.3f P=%.3g\n",
            sum(m$dis == "OA"), tp_oa$estimate, tp_oa$p.value))
cat(sprintf("[A2] Medial vs Lateral 仅Normal供体(n=%d对): %+.3f P=%.3g\n",
            sum(m$dis == "Normal"), tp_no$estimate, tp_no$p.value))
addrow(Dataset = "GSE51588", Tissue = "Subchondral bone",
       Comparison = "Load-bearing Medial vs non-load Lateral (all donors)",
       Probe = pz, n_Ctrl = nrow(m), n_Case = nrow(m),
       Diff = tp_all$estimate, CI_low = tp_all$conf.int[1], CI_high = tp_all$conf.int[2],
       P = tp_all$p.value, AUC = NA, Method = "Paired t-test", Direction = "Medial - Lateral")
addrow(Dataset = "GSE51588", Tissue = "Subchondral bone",
       Comparison = "Medial vs Lateral (OA donors only)",
       Probe = pz, n_Ctrl = sum(m$dis == "OA"), n_Case = sum(m$dis == "OA"),
       Diff = tp_oa$estimate, CI_low = tp_oa$conf.int[1], CI_high = tp_oa$conf.int[2],
       P = tp_oa$p.value, AUC = NA, Method = "Paired t-test", Direction = "Medial - Lateral")
addrow(Dataset = "GSE51588", Tissue = "Subchondral bone",
       Comparison = "Medial vs Lateral (Normal donors only)",
       Probe = pz, n_Ctrl = sum(m$dis == "Normal"), n_Case = sum(m$dis == "Normal"),
       Diff = tp_no$estimate, CI_low = tp_no$conf.int[1], CI_high = tp_no$conf.int[2],
       P = tp_no$p.value, AUC = NA, Method = "Paired t-test", Direction = "Medial - Lateral")

## ---- A3 疾病 x 应力交互: delta 在 OA vs Normal 间 ----
t3 <- t.test(delta ~ dis, data = m)
d3i <- mean(m$delta[m$dis == "OA"]) - mean(m$delta[m$dis == "Normal"])
cat(sprintf("[A3] 交互(delta=Medial-Lateral): OA %+.3f vs Normal %+.3f | P=%.3g\n",
            mean(m$delta[m$dis == "OA"]), mean(m$delta[m$dis == "Normal"]), t3$p.value))
addrow(Dataset = "GSE51588", Tissue = "Subchondral bone",
       Comparison = "Interaction: delta(Medial-Lateral) OA vs Normal",
       Probe = pz, n_Ctrl = sum(m$dis == "Normal"), n_Case = sum(m$dis == "OA"),
       Diff = d3i, CI_low = t3$conf.int[1], CI_high = t3$conf.int[2],
       P = t3$p.value, AUC = NA, Method = "Two-sample t-test on paired deltas",
       Direction = "OA - Normal")

## ---- A4 年龄相关 ----
ok <- !is.na(age)
ct_all <- cor.test(v[ok], age[ok], method = "spearman")
ct_oa  <- cor.test(v[ok & dis == "OA"], age[ok & dis == "OA"], method = "spearman")
cat(sprintf("[A4] 年龄相关(Spearman): 全部 rho=%.2f P=%.3g | 仅OA rho=%.2f P=%.3g\n",
            ct_all$estimate, ct_all$p.value, ct_oa$estimate, ct_oa$p.value))

## ---- 图1: 四分组箱线图 ----
df4 <- data.frame(group = g4, expr = v)
lv4 <- c("Normal\nLateral", "Normal\nMedial", "OA\nLateral", "OA\nMedial")
df4$group <- factor(df4$group, levels = lv4)
p1 <- ggboxplot(df4, x = "group", y = "expr", fill = "group", add = "jitter",
                add.params = list(size = 1.6, alpha = 0.7), palette = "npg",
                xlab = NULL, ylab = "PIEZO1 expression (normalized)",
                title = "PIEZO1 in subchondral bone - GSE51588",
                legend = "none") +
  stat_compare_means(comparisons = list(
    c("Normal\nLateral", "OA\nLateral"), c("Normal\nMedial", "OA\nMedial"),
    c("Normal\nLateral", "Normal\nMedial"), c("OA\nLateral", "OA\nMedial")),
    method = "t.test", size = 2.6) +
  rotate_x_text(0) +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(figs, "fig01_GSE51588_fourgroups.png"), p1,
       width = 6.2, height = 5.4, dpi = 300)

## ---- 图2: 配对折线图 (同一供体 内侧 vs 外侧) ----
dfp <- w; dfp$reg <- factor(dfp$reg, levels = c("Lateral", "Medial"))
dfp$dis <- factor(dfp$dis, levels = c("Normal", "OA"))
p2 <- ggplot(dfp, aes(x = reg, y = expr, group = donor)) +
  geom_line(color = "grey55", linewidth = 0.55, alpha = 0.85) +
  geom_point(aes(color = dis), size = 2.1) +
  facet_wrap(~ dis, ncol = 2) +
  scale_color_manual(values = c(Normal = BLUE, OA = RED)) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", size = 3.2,
               shape = 18, color = "black") +
  stat_summary(aes(group = 1), fun = mean, geom = "line",
               linewidth = 1.2, color = "black", linetype = "dashed") +
  labs(title = "PIEZO1: load-bearing (Medial) vs non-load (Lateral) within the same donor",
       subtitle = "Each line = one donor; diamond = group mean (GSE51588)",
       x = NULL, y = "PIEZO1 expression (normalized)", color = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 11.5, face = "bold"),
        strip.text = element_text(face = "bold"))
ggsave(file.path(figs, "fig02_GSE51588_paired_lines.png"), p2,
       width = 7.6, height = 4.6, dpi = 300)

## ---- 图3: 机械应力敏感性 delta 箱线图 ----
p3 <- ggboxplot(m, x = "dis", y = "delta", fill = "dis", add = "jitter",
                add.params = list(size = 2, alpha = 0.7), palette = c(BLUE, RED),
                xlab = NULL, ylab = "PIEZO1 delta = Medial - Lateral",
                title = "Mechanical-stress responsiveness of PIEZO1 (GSE51588)",
                legend = "none") +
  stat_compare_means(method = "t.test", size = 3.2, label.x.npc = "center") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(figs, "fig03_GSE51588_mechanical_delta.png"), p3,
       width = 4.8, height = 5.0, dpi = 300)

## ---- 图4: 年龄相关散点 ----
dfa <- data.frame(age = age, expr = v, dis = dis)[ok, ]
p4 <- ggscatter(dfa, x = "age", y = "expr", color = "dis",
               add = "reg.line", conf.int = TRUE,
               palette = c(BLUE, RED),
               xlab = "Age (years)", ylab = "PIEZO1 expression",
               title = "PIEZO1 vs age - GSE51588") +
  stat_cor(aes(color = dis), method = "spearman", size = 2.8) +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(figs, "fig04_GSE51588_age_correlation.png"), p4,
       width = 6.0, height = 4.6, dpi = 300)

## ---- 图5: ROC (OA vs Normal) ----
png(file.path(figs, "fig05_GSE51588_ROC.png"), width = 1600, height = 1600, res = 300)
par(mar = c(4.2, 4.2, 3, 1))
plot(roc1, col = RED, lwd = 2.5, legacy.axes = TRUE,
     main = sprintf("ROC of PIEZO1: OA vs Normal (GSE51588 subchondral bone)\nAUC = %.3f", auc1),
     xlab = "1 - Specificity", ylab = "Sensitivity")
abline(0, 1, lty = 2, col = "grey60")
dev.off()

# ============================================================
# PART B: GSE82107 滑膜 (多探针一致性)
# ============================================================
cat("\n========== B. GSE82107 滑膜 ==========\n")
g2 <- getGEO(filename = file.path(d, "data", "GSE82107_series_matrix.txt.gz"),
             getGPL = FALSE)
e2 <- as.matrix(exprs(g2)); pd2 <- pData(g2)
if (max(e2, na.rm = TRUE) > 100) e2 <- log2(e2 + 1)
grp_raw <- NULL
for (cc in colnames(pd2)) {
  vv <- as.character(pd2[[cc]])
  if (any(grepl("osteoarthrit", vv, ignore.case = TRUE))) { grp_raw <- vv; break }
}
vv2 <- toupper(trimws(sub(".*:", "", grp_raw)))
gr2 <- ifelse(vv2 == "OA" | grepl("OSTEOARTH", vv2), "OA", "Healthy")
gr2 <- factor(gr2, levels = c("Healthy", "OA"))
probes82107 <- c("202771_at", "1566110_at", "1566111_at")
cat(sprintf("分组: %s\n", paste(names(table(gr2)), table(gr2), sep = "=", collapse = ", ")))
stat82107 <- NULL
for (pb in probes82107) {
  if (!pb %in% rownames(e2)) next
  x <- as.numeric(e2[pb, ])
  dd <- mean(x[gr2 == "OA"]) - mean(x[gr2 == "Healthy"])   # OA - Healthy
  tt <- t.test(x[gr2 == "OA"], x[gr2 == "Healthy"])
  rr <- roc(response = as.integer(gr2 == "OA"), predictor = x,
            quiet = TRUE, direction = "auto")
  stat82107 <- rbind(stat82107, data.frame(probe = pb, diff = dd,
                                           P = tt$p.value, AUC = as.numeric(auc(rr))))
  cat(sprintf("[B] %s: %+.3f P=%.3g AUC=%.3f\n", pb, dd,
              tt$p.value, as.numeric(auc(rr))))
  if (pb == "202771_at") {
    addrow(Dataset = "GSE82107", Tissue = "Synovium", Comparison = "OA vs Healthy",
           Probe = pb, n_Ctrl = sum(gr2 == "Healthy"), n_Case = sum(gr2 == "OA"),
           Diff = dd, CI_low = tt$conf.int[1], CI_high = tt$conf.int[2],
           P = tt$p.value, AUC = as.numeric(auc(rr)), Method = "Welch t-test",
           Direction = "OA - Healthy")
    roc2 <- rr
  }
}
## ---- 图6: 三探针 facet 箱线图 ----
dfp2 <- do.call(rbind, lapply(probes82107[probes82107 %in% rownames(e2)], function(pb) {
  data.frame(probe = pb, group = gr2, expr = as.numeric(e2[pb, ]))
}))
p6 <- ggboxplot(dfp2, x = "group", y = "expr", fill = "group", add = "jitter",
                add.params = list(size = 1.4, alpha = 0.7), palette = c(BLUE, RED),
                facet.by = "probe", scales = "free_y", xlab = NULL,
                ylab = "PIEZO1 expression (log2)",
                title = "PIEZO1 probes in OA synovium - GSE82107") +
  stat_compare_means(method = "t.test", size = 3, label.y.npc = 0.95) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(face = "bold"))
ggsave(file.path(figs, "fig06_GSE82107_probes.png"), p6,
       width = 8.2, height = 4.8, dpi = 300)
## ---- 图7: ROC (主探针 202771_at) ----
png(file.path(figs, "fig07_GSE82107_ROC.png"), width = 1600, height = 1600, res = 300)
par(mar = c(4.2, 4.2, 3, 1))
plot(roc2, col = RED, lwd = 2.5, legacy.axes = TRUE,
     main = sprintf("ROC of PIEZO1 (202771_at): OA vs Healthy\nGSE82107 synovium   AUC = %.3f",
                    as.numeric(auc(roc2))),
     xlab = "1 - Specificity", ylab = "Sensitivity")
abline(0, 1, lty = 2, col = "grey60")
dev.off()

# ============================================================
# PART C: GSE55235 滑膜 (统一风格重画)
# ============================================================
cat("\n========== C. GSE55235 滑膜 ==========\n")
g3 <- getGEO(filename = data_path("GSE55235_series_matrix.txt.gz"),
             getGPL = FALSE)
e3 <- as.matrix(exprs(g3)); pd3 <- pData(g3)
if (max(e3, na.rm = TRUE) > 100) e3 <- log2(e3 + 1)
chr_cols <- grep("characteristics", colnames(pd3), value = TRUE)
grp_raw3 <- NULL
for (cc in chr_cols) {
  vv <- as.character(pd3[[cc]])
  key <- trimws(sub(":.*", "", vv))
  if (length(unique(key)) == 1 && grepl("disease|status|condition|state", key[1], ignore.case = TRUE)) {
    grp_raw3 <- trimws(sub(".*:", "", vv)); names(grp_raw3) <- rownames(pd3); break
  }
}
is_case <- grepl("osteoarthrit", grp_raw3, ignore.case = TRUE)
is_norm <- grepl("healthy|normal|control", grp_raw3, ignore.case = TRUE)
keep <- is_case | is_norm
gr3 <- factor(ifelse(is_case, "OA", "Healthy")[keep], levels = c("Healthy", "OA"))
names(gr3) <- names(grp_raw3)[keep]
e3 <- e3[, names(gr3)]
pb3 <- "202771_at"
x3 <- as.numeric(e3[pb3, ])
d3 <- mean(x3[gr3 == "OA"]) - mean(x3[gr3 == "Healthy"])   # OA - Healthy
t3b <- t.test(x3[gr3 == "OA"], x3[gr3 == "Healthy"])
roc3 <- roc(response = as.integer(gr3 == "OA"), predictor = x3,
            quiet = TRUE, direction = "auto")
auc3 <- as.numeric(auc(roc3))
cat(sprintf("[C] %s: %+.3f (95%%CI %+.3f~%+.3f) P=%.3g AUC=%.3f\n", pb3,
            d3, t3b$conf.int[1], t3b$conf.int[2], t3b$p.value, auc3))
addrow(Dataset = "GSE55235", Tissue = "Synovium", Comparison = "OA vs Healthy",
       Probe = pb3, n_Ctrl = sum(gr3 == "Healthy"), n_Case = sum(gr3 == "OA"),
       Diff = d3, CI_low = t3b$conf.int[1], CI_high = t3b$conf.int[2],
       P = t3b$p.value, AUC = auc3, Method = "Welch t-test", Direction = "OA - Healthy")
## ---- 图8: 箱线图 ----
df3 <- data.frame(group = gr3, expr = x3)
p8 <- ggboxplot(df3, x = "group", y = "expr", fill = "group", add = "jitter",
                add.params = list(size = 1.8, alpha = 0.7), palette = c(BLUE, RED),
                xlab = NULL, ylab = "PIEZO1 expression (log2)",
                title = "PIEZO1 in OA synovium - GSE55235 (202771_at)",
                legend = "none") +
  stat_compare_means(method = "t.test", size = 3.4) +
  theme(plot.title = element_text(size = 12, face = "bold"))
ggsave(file.path(figs, "fig08_GSE55235_boxplot.png"), p8,
       width = 4.8, height = 5.0, dpi = 300)
## ---- 图9: ROC ----
png(file.path(figs, "fig09_GSE55235_ROC.png"), width = 1600, height = 1600, res = 300)
par(mar = c(4.2, 4.2, 3, 1))
plot(roc3, col = RED, lwd = 2.5, legacy.axes = TRUE,
     main = sprintf("ROC of PIEZO1 (202771_at): OA vs Healthy\nGSE55235 synovium   AUC = %.3f", auc3),
     xlab = "1 - Specificity", ylab = "Sensitivity")
abline(0, 1, lty = 2, col = "grey60")
dev.off()

# ============================================================
# PART D: 跨数据集森林图 + 汇总表
# ============================================================
cat("\n========== D. 汇总 ==========\n")
ff <- rbindlist(sumrows, fill = TRUE)
setorder(ff, Dataset, Comparison)
fwrite(ff, file.path(res, "PIEZO1_expression_diagnosis_summary.csv"))
print(ff[, .(Dataset, Comparison, Diff = round(Diff, 3),
            CI_low = round(CI_low, 3), CI_high = round(CI_high, 3),
            P = signif(P, 3), AUC = round(AUC, 3))])
ff[, label := paste0(Dataset, " | ", Comparison)]
ff[, sig := ifelse(P < 0.05, "P < 0.05", "NS")]
ff[, ptxt := sprintf("P=%.3g", P)]
ff[, lab := factor(seq_len(nrow(ff)), levels = rev(seq_len(nrow(ff))))]
## ---- 图10: 森林图 ----
p10 <- ggplot(ff, aes(x = Diff, y = lab)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.22,
                 linewidth = 0.7, color = "grey35") +
  geom_point(aes(color = Dataset, shape = sig), size = 3.4) +
  geom_text(aes(x = max(CI_high, na.rm = TRUE), label = ptxt), hjust = -0.15,
            size = 3, color = "grey20") +
  scale_shape_manual(values = c(`P < 0.05` = 16, NS = 1),
                     name = "Significance") +
  scale_color_brewer(palette = "Set2", name = "Dataset") +
  labs(title = "PIEZO1 in osteoarthritis: summary of comparisons",
       subtitle = "Mean difference with 95% CI; filled = P < 0.05",
       x = "Mean difference (case - control)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(lab ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(lab = function(x) as.character(ff$label[x]))) +
  theme(strip.text.y.left = element_text(angle = 0, size = 8.6, hjust = 1))
ggsave(file.path(figs, "fig10_forest_all.png"), p10,
       width = 8.6, height = 5.6, dpi = 300)
cat("\n全部完成: 10 张图 + 1 个汇总表\n")
cat("图: ", paste(sprintf("fig%02d", 1:10), collapse = ", "), "\n")
