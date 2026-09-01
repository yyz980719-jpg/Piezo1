# ==============================================================================
# 40_P0_donorlevel_inference.R
# P0-2: 修复 GSE104782 的假重复问题（审稿人 R1 第 2 点 / 交叉综合共识风险 #2）
#
# 问题：稿件头条数字 Kruskal-Wallis P = 3.9e-6 与 Spearman rho = -0.142, P = 5.3e-8
#       均为细胞水平检验，但 1,464 个细胞嵌套于 10 个供体，细胞并非独立生物学重复；
#       且 rho = -0.142 仅解释约 2% 方差，P 值主要由大 N 驱动。
#
# 三条互为补充的供体级推断路径：
#   (A) pseudobulk：按 patient x subtype 聚合计数 -> CPM -> limma（供体为区组）
#   (B) 细胞水平 limma + duplicateCorrelation(block = patient)（供体为随机效应）
#   (C) 供体内配对 + 跨供体检验（每供体算一个统计量，再对 10 个供体做推断）
#
# 输出: results/40_*.csv, figs/fig60_*
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(limma); library(edgeR)
  library(data.table); library(ggplot2); library(patchwork); library(dplyr)
})

set.seed(20260830)
source("scripts/_bootstrap.R")
RES <- PIEZO1_RESULTS; FIG <- PIEZO1_FIGURES
dir.create(RES, showWarnings = FALSE); dir.create(FIG, showWarnings = FALSE)

theme_set(theme_bw(base_size = 11) +
            theme(panel.grid = element_blank(),
                  plot.title = element_text(face = "bold", size = 11)))

sink(file.path(RES, "40_P0_donorlevel_log.txt"), split = TRUE)
cat("================================================================\n")
cat(" P0-2  GSE104782 供体级推断（修复假重复）\n")
cat(" 运行时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(" R:", R.version.string, "| limma:", as.character(packageVersion("limma")),
    "| Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("================================================================\n\n")

## ---------- 0. 读入 ----------
seu <- readRDS("data/scrna/GSE104782_seurat.rds")
md  <- seu@meta.data
md$subtype <- as.character(md$subtype)
md$patient <- as.character(md$patient)

cnt <- as.matrix(GetAssayData(seu, layer = "counts"))
nrm <- as.matrix(GetAssayData(seu, layer = "data"))

cat(sprintf("细胞数 %d | 基因数 %d | 供体 %d | 亚型 %d\n",
            ncol(cnt), nrow(cnt), length(unique(md$patient)), length(unique(md$subtype))))
stopifnot("PIEZO1" %in% rownames(cnt))

tab_ps <- table(md$patient, md$subtype)
cat("\n--- patient x subtype 细胞数 ---\n"); print(tab_ps)
cat("\n每亚型总细胞数:\n"); print(colSums(tab_ps))
cat("\n每亚型有细胞的供体数:\n"); print(colSums(tab_ps > 0))
cat("\n有 >=5 个细胞的 donor x subtype 组合数（pseudobulk 可用）:\n")
print(colSums(tab_ps >= 5))
tab_out <- as.data.frame.matrix(tab_ps)
fwrite(data.frame(patient = rownames(tab_out), tab_out, row.names = NULL),
       file.path(RES, "40_P0_cells_per_donor_subtype.csv"))

## ---------- 1. 细胞水平：复现稿件现有数字（对照基线） ----------
pz_cell <- as.numeric(nrm["PIEZO1", ])
kw_cell <- kruskal.test(pz_cell ~ factor(md$subtype))

## HGNC/NCBI canonical symbol: HAPLN1. CRTL1 is a historical alias of HAPLN1
## and must not be counted as an additional marker.
homeo <- c("COL2A1","SOX9","ACAN","CHAD","HAPLN1")
hyp   <- c("COL10A1","RUNX2","IBSP","ALPL","MMP13","SPP1","POSTN")
homeo <- intersect(homeo, rownames(nrm)); hyp <- intersect(hyp, rownames(nrm))
sc <- function(gs) as.numeric(scale(colMeans(nrm[gs, , drop = FALSE])))
axis <- sc(hyp) - sc(homeo)                       # 正 = 肥大极（与 38_diff_axis.R 一致）
cat(sprintf("\n分化轴 marker: homeo = %d, hyp = %d\n", length(homeo), length(hyp)))

ct_cell  <- cor.test(pz_cell, axis, method = "spearman")
q_all    <- cut(axis, breaks = quantile(axis, probs = seq(0, 1, 0.25)),
                include.lowest = TRUE, labels = c("Q1_homeo","Q2","Q3","Q4_hyp"))
prof_cell <- tapply(pz_cell, q_all, mean)

cat("\n########## 【基线】细胞水平（稿件现状，存在假重复）##########\n")
cat(sprintf("Kruskal-Wallis (7 亚型):         P = %.3g\n", kw_cell$p.value))
cat(sprintf("PIEZO1 vs 分化轴 Spearman:       rho = %.3f, P = %.3g\n",
            ct_cell$estimate, ct_cell$p.value))
cat(sprintf("四分位均值 Q1(稳态) -> Q4(肥大): %.3f -> %.3f\n",
            prof_cell[["Q1_homeo"]], prof_cell[["Q4_hyp"]]))
cat(sprintf("效应量: rho = %.3f ≈ %.1f%% 方差解释率 —— P 值由 n = %d 个细胞驱动\n",
            ct_cell$estimate, 100 * ct_cell$estimate^2, length(pz_cell)))

## ---------- 2. 路径 A：pseudobulk（patient x subtype）+ limma ----------
cat("\n########## 【路径 A】pseudobulk: patient x subtype ##########\n")

grp     <- interaction(md$patient, md$subtype, sep = "__", drop = TRUE)
pb_cnt  <- sapply(levels(grp), function(g) Matrix::rowSums(cnt[, grp == g, drop = FALSE]))
pm      <- do.call(rbind, strsplit(levels(grp), "__", fixed = TRUE))
pb_meta <- data.frame(patient = pm[, 1], subtype = pm[, 2], stringsAsFactors = FALSE)
pb_meta$n_cells <- as.integer(tab_ps[cbind(pb_meta$patient, pb_meta$subtype)])

MIN_CELLS <- 5
keep_pb   <- pb_meta$n_cells >= MIN_CELLS
cat(sprintf("pseudobulk 组合: %d 个, 保留 n_cells >= %d 的 %d 个\n",
            length(keep_pb), MIN_CELLS, sum(keep_pb)))
pb_cnt  <- pb_cnt[, keep_pb, drop = FALSE]
pb_meta <- pb_meta[keep_pb, ]
cat("保留组合 — 供体分布:\n"); print(table(pb_meta$patient))
cat("保留组合 — 亚型分布:\n"); print(table(pb_meta$subtype))

lib <- colSums(pb_cnt)
cat(sprintf("pseudobulk 文库大小: 中位 %.0f (范围 %.0f – %.0f)\n",
            median(lib), min(lib), max(lib)))

cpm_mat <- edgeR::cpm(pb_cnt, log = TRUE, prior.count = 2)

## HomC 设为参照水平；~ patient + subtype（含截距，保证满秩）
sub_levs <- c("HomC", setdiff(sort(unique(pb_meta$subtype)), "HomC"))
pb_meta$patient <- factor(pb_meta$patient)
pb_meta$subtype <- factor(pb_meta$subtype, levels = sub_levs)
design_A <- model.matrix(~ patient + subtype, data = pb_meta)
cat(sprintf("\n设计矩阵: %d 样本 x %d 系数 (满秩 = %d)\n",
            nrow(design_A), ncol(design_A), qr(design_A)$rank))

fitA <- lmFit(cpm_mat, design_A)
fitA <- eBayes(fitA, trend = TRUE)

sub_coefs <- grep("^subtype", colnames(design_A), value = TRUE)
cat("亚型系数（相对 HomC 参照）:", paste(sub_coefs, collapse = ", "), "\n")

## omnibus F 检验（供体校正后，全部亚型效应）
omnA <- topTable(fitA, coef = sub_coefs, number = Inf)
cat("\n--- PIEZO1 全亚型 omnibus F 检验（供体区组校正）---\n")
cat(sprintf("F = %.3f | P = %.3g | FDR = %.3g\n",
            omnA["PIEZO1","F"],
            omnA["PIEZO1","P.Value"], omnA["PIEZO1","adj.P.Val"]))
piezo1_idx <- match("PIEZO1", rownames(fitA$coefficients))
omnibus_df1 <- length(sub_coefs)
omnibus_df2 <- fitA$df.total[piezo1_idx]
cat(sprintf("Moderated omnibus degrees of freedom: df1 = %d | df2 = %.6f\n",
            omnibus_df1, omnibus_df2))

## 关键对比
cmA <- makeContrasts(
  HomC_minus_HTC    = -subtypeHTC,
  HomC_minus_FC     = -subtypeFC,
  HomC_minus_preHTC = -subtypepreHTC,
  HomC_minus_ProC   = -subtypeProC,
  RegC_minus_HTC    =  subtypeRegC - subtypeHTC,
  RegC_minus_FC     =  subtypeRegC - subtypeFC,
  levels = design_A)
cmA <- cmA[, apply(cmA, 2, function(z) all(names(z[abs(z) > 0]) %in% colnames(design_A))), drop = FALSE]

resA <- do.call(rbind, lapply(colnames(cmA), function(cc) {
  ttc <- topTable(contrasts.fit(fitA, cmA[, cc, drop = FALSE]) |> eBayes(),
                  coef = 1, number = Inf, confint = TRUE)
  data.frame(contrast = cc, logFC = round(ttc["PIEZO1","logFC"], 4),
             CI.L = round(ttc["PIEZO1","CI.L"], 4), CI.R = round(ttc["PIEZO1","CI.R"], 4),
             P = signif(ttc["PIEZO1","P.Value"], 4), FDR = signif(ttc["PIEZO1","adj.P.Val"], 4))
}))
cat("\n--- 关键两两对比（pseudobulk，供体区组校正）---\n")
print(resA, row.names = FALSE)
fwrite(resA, file.path(RES, "40_P0_pseudobulk_piezo1_contrasts.csv"))
omnA_all <- data.frame(gene = rownames(omnA), omnA[, c("F","P.Value","adj.P.Val")],
                       row.names = NULL)
fwrite(omnA_all, file.path(RES, "40_P0_pseudobulk_omnibus_all_genes.csv"))
fwrite(data.frame(gene = "PIEZO1", F = omnA["PIEZO1","F"],
                  df1 = omnibus_df1, df2 = omnibus_df2,
                  P = omnA["PIEZO1","P.Value"], FDR = omnA["PIEZO1","adj.P.Val"]),
       file.path(RES, "40_P0_pseudobulk_piezo1_omnibus.csv"))

## ---------- 3. 路径 B：细胞水平 limma + duplicateCorrelation ----------
cat("\n########## 【路径 B】细胞水平 limma + duplicateCorrelation(block = patient) ##########\n")

## 全转录组 duplicateCorrelation 对本任务开销很大；主推断已由 donor-blocked
## pseudobulk 和下方供体级配对检验覆盖。仅在显式需要时开启。
RUN_DUPCOR <- FALSE

if (RUN_DUPCOR) {

mdB <- data.frame(patient = factor(md$patient), subtype = factor(md$subtype))
design_B <- model.matrix(~ 0 + subtype, data = mdB)
cat(sprintf("细胞水平设计: %d 细胞 x %d 系数, 区组 = %d 供体\n",
            nrow(design_B), ncol(design_B), nlevels(mdB$patient)))

dupcor <- duplicateCorrelation(nrm, design_B, block = mdB$patient)
cat(sprintf("供体内相关 (consensus correlation) = %.4f\n", dupcor$consensus.correlation))

fitB <- lmFit(nrm, design_B, block = mdB$patient,
              correlation = dupcor$consensus.correlation)
cmB <- makeContrasts(
  HomC_minus_HTC    = subtypeHomC - subtypeHTC,
  HomC_minus_FC     = subtypeHomC - subtypeFC,
  HomC_minus_preHTC = subtypeHomC - subtypepreHTC,
  RegC_minus_HTC    = subtypeRegC - subtypeHTC,
  HomC_minus_ProC   = subtypeHomC - subtypeProC,
  levels = design_B)
cmB <- cmB[, apply(cmB, 2, function(z) all(names(z[abs(z) > 0]) %in% colnames(design_B))), drop = FALSE]

fitBc <- contrasts.fit(fitB, cmB) |> eBayes()
resB <- do.call(rbind, lapply(colnames(cmB), function(cc) {
  ttc <- topTable(fitBc, coef = cc, number = Inf, confint = TRUE)
  data.frame(contrast = cc, logFC = round(ttc["PIEZO1","logFC"], 4),
             CI.L = round(ttc["PIEZO1","CI.L"], 4), CI.R = round(ttc["PIEZO1","CI.R"], 4),
             P = signif(ttc["PIEZO1","P.Value"], 4), FDR = signif(ttc["PIEZO1","adj.P.Val"], 4))
}))
cat("\n--- 细胞水平 + 供体随机效应校正（供体内相关已纳入）---\n")
print(resB, row.names = FALSE)
fwrite(resB, file.path(RES, "40_P0_duplicatecorrelation_piezo1.csv"))
fwrite(data.frame(consensus_correlation = dupcor$consensus.correlation),
       file.path(RES, "40_P0_duplicatecorrelation_value.csv"))
} else {
  cat("跳过：使用 donor-blocked pseudobulk + donor-level paired inference 作为预设替代。\n")
}

## ---------- 4. 路径 C：供体内配对 + 跨供体推断（分化轴） ----------
cat("\n########## 【路径 C】供体内配对 + 跨供体推断（分化轴趋势）##########\n")

donors <- sort(unique(md$patient))
donor_tab <- do.call(rbind, lapply(donors, function(d) {
  i <- md$patient == d
  x <- pz_cell[i]; a <- axis[i]; q <- q_all[i]
  ct <- suppressWarnings(cor.test(x, a, method = "spearman"))
  data.frame(patient = d, n_cells = sum(i),
             rho = unname(ct$estimate), P_cell = ct$p.value,
             mean_Q1_homeo = if (sum(q == "Q1_homeo") > 1) mean(x[q == "Q1_homeo"]) else NA,
             mean_Q4_hyp   = if (sum(q == "Q4_hyp")   > 1) mean(x[q == "Q4_hyp"])   else NA,
             n_Q1 = sum(q == "Q1_homeo"), n_Q4 = sum(q == "Q4_hyp"),
             stringsAsFactors = FALSE)
}))
donor_tab$delta_Q1_minus_Q4 <- donor_tab$mean_Q1_homeo - donor_tab$mean_Q4_hyp
donor_tab$z_rho <- atanh(donor_tab$rho)

cat("\n--- 每位供体的分化轴关联 ---\n")
print(data.frame(patient = donor_tab$patient, n_cells = donor_tab$n_cells,
                 rho = round(donor_tab$rho, 3), P_cell = signif(donor_tab$P_cell, 3),
                 Q1 = round(donor_tab$mean_Q1_homeo, 3),
                 Q4 = round(donor_tab$mean_Q4_hyp, 3),
                 delta = round(donor_tab$delta_Q1_minus_Q4, 3)), row.names = FALSE)

t_rho <- t.test(donor_tab$rho)
t_z   <- t.test(donor_tab$z_rho)
w_rho <- suppressWarnings(wilcox.test(donor_tab$rho, exact = FALSE))
n_neg <- sum(donor_tab$rho < 0, na.rm = TRUE)
n_tot <- sum(!is.na(donor_tab$rho))
sign_p <- binom.test(n_neg, n_tot, p = 0.5)$p.value

ok  <- !is.na(donor_tab$delta_Q1_minus_Q4)
t_q <- t.test(donor_tab$mean_Q1_homeo[ok], donor_tab$mean_Q4_hyp[ok], paired = TRUE)
w_q <- suppressWarnings(wilcox.test(donor_tab$mean_Q1_homeo[ok], donor_tab$mean_Q4_hyp[ok],
                                    paired = TRUE, exact = FALSE))

cat("\n--- 供体级推断（n = 10 供体，正确推断单元）---\n")
cat(sprintf("供体水平 rho: 均值 = %.3f, 中位 = %.3f, 范围 [%.3f, %.3f]\n",
            mean(donor_tab$rho), median(donor_tab$rho),
            min(donor_tab$rho), max(donor_tab$rho)))
cat(sprintf("单样本 t 检验 (rho):       t = %.3f, df = %g, P = %.4g, 95%%CI [%.3f, %.3f]\n",
            t_rho$statistic, t_rho$parameter, t_rho$p.value,
            t_rho$conf.int[1], t_rho$conf.int[2]))
cat(sprintf("Fisher-z 单样本 t 检验:    P = %.4g\n", t_z$p.value))
cat(sprintf("Wilcoxon 符号秩检验 (rho): P = %.4g\n", w_rho$p.value))
cat(sprintf("符号检验: %d/%d 供体 rho < 0, P = %.4g\n", n_neg, n_tot, sign_p))
cat(sprintf("\nQ1(稳态) vs Q4(肥大) 配对 t 检验: 均值差 = %.4f, P = %.4g, 95%%CI [%.4f, %.4f]\n",
            t_q$estimate, t_q$p.value, t_q$conf.int[1], t_q$conf.int[2]))
cat(sprintf("Q1 vs Q4 配对 Wilcoxon:            P = %.4g\n", w_q$p.value))

fwrite(donor_tab, file.path(RES, "40_P0_donorlevel_axis_perdonor.csv"))

axis_tests <- data.frame(
  test = c("one-sample t on donor rho", "Fisher-z one-sample t",
           "Wilcoxon signed-rank (donor rho)", "sign test (rho<0)",
           "paired t Q1 vs Q4", "Wilcoxon paired Q1 vs Q4",
           "cell-level Spearman (for comparison)"),
  statistic = c(unname(t_rho$statistic), unname(t_z$statistic), unname(w_rho$statistic),
                n_neg, unname(t_q$statistic), unname(w_q$statistic), unname(ct_cell$estimate)),
  n_unit = c(n_tot, n_tot, n_tot, n_tot, sum(ok), sum(ok), length(pz_cell)),
  unit = c("donors","donors","donors","donors","donors","donors","cells"),
  P = signif(c(t_rho$p.value, t_z$p.value, w_rho$p.value, sign_p,
               t_q$p.value, w_q$p.value, ct_cell$p.value), 4),
  effect = c(sprintf("mean rho = %.3f", mean(donor_tab$rho)),
             sprintf("mean z = %.3f", mean(donor_tab$z_rho)),
             sprintf("median rho = %.3f", median(donor_tab$rho)),
             sprintf("%d/%d negative", n_neg, n_tot),
             sprintf("mean diff = %.4f", t_q$estimate),
             sprintf("median diff = %.4f", median(donor_tab$delta_Q1_minus_Q4[ok])),
             sprintf("rho = %.3f (%.1f%% var)", ct_cell$estimate, 100*ct_cell$estimate^2)))
fwrite(axis_tests, file.path(RES, "40_P0_axis_donorlevel_tests.csv"))

## ---------- 5. 路径 C 扩展：亚型差异的供体级（配对）检验 ----------
cat("\n########## 【路径 C 扩展】亚型差异的供体级（配对）检验 ##########\n")

pz_mat <- tapply(pz_cell, list(md$patient, md$subtype), mean)   # 供体 x 亚型
cat("\n--- 供体 x 亚型 PIEZO1 平均表达 ---\n"); print(round(pz_mat, 3))
fwrite(data.frame(patient = rownames(pz_mat), round(pz_mat, 4), row.names = NULL),
       file.path(RES, "40_P0_donor_x_subtype_piezo1.csv"))

n_have <- colSums(!is.na(pz_mat))
cat("\n每亚型有数据的供体数:\n"); print(n_have)
subs_full <- names(n_have)[n_have >= 8]
cat(sprintf("可用于配对检验的亚型 (>=8 供体): %s\n", paste(subs_full, collapse = ", ")))

pzm <- pz_mat[, subs_full, drop = FALSE]
cc  <- complete.cases(pzm)
cat(sprintf("完整供体（所选亚型均有值）: %d 个\n", sum(cc)))

if (sum(cc) >= 5 && ncol(pzm) >= 3) {
  fr <- friedman.test(as.matrix(pzm[cc, ]))
  W  <- unname(fr$statistic) / (sum(cc) * (ncol(pzm) - 1))
  cat(sprintf("Friedman 检验 (%d 供体 x %d 亚型): chi2 = %.3f, df = %g, P = %.4g\n",
              sum(cc), ncol(pzm), fr$statistic, fr$parameter, fr$p.value))
  cat(sprintf("Kendall's W (效应量) = %.3f\n", W))
  fr_out <- data.frame(test = "Friedman", statistic = unname(fr$statistic),
                       df = unname(fr$parameter), P = signif(fr$p.value, 4),
                       n_donors = sum(cc), n_subtypes = ncol(pzm), kendall_W = round(W, 3))
} else {
  fr_out <- data.frame(test = "Friedman", statistic = NA, df = NA, P = NA,
                       n_donors = sum(cc), n_subtypes = ncol(pzm), kendall_W = NA)
}
fwrite(fr_out, file.path(RES, "40_P0_donorlevel_subtype_omnibus.csv"))

pair_list <- list(c("HomC","HTC"), c("HomC","FC"), c("RegC","HTC"),
                  c("HomC","preHTC"), c("RegC","FC"), c("HomC","RegC"),
                  c("HomC","EC"), c("HomC","ProC"))
pair_out <- do.call(rbind, lapply(pair_list, function(pp) {
  if (!all(pp %in% colnames(pz_mat))) return(NULL)
  a <- pz_mat[, pp[1]]; b <- pz_mat[, pp[2]]
  k <- !is.na(a) & !is.na(b)
  if (sum(k) < 4) return(NULL)
  tt <- t.test(a[k], b[k], paired = TRUE)
  ww <- suppressWarnings(wilcox.test(a[k], b[k], paired = TRUE, exact = FALSE))
  d  <- mean(a[k] - b[k]); s <- sd(a[k] - b[k])
  data.frame(contrast = paste(pp, collapse = " vs "), n_donors = sum(k),
             mean_diff = round(d, 4), cohens_dz = round(d / s, 3),
             CI_low = round(tt$conf.int[1], 4), CI_high = round(tt$conf.int[2], 4),
             P_paired_t = signif(tt$p.value, 4), P_wilcoxon = signif(ww$p.value, 4))
}))
cat("\n--- 供体内配对检验（亚型均值，n = 供体数）---\n")
print(pair_out, row.names = FALSE)
fwrite(pair_out, file.path(RES, "40_P0_donorlevel_subtype_paired.csv"))

## ---------- 6. 汇总对照表 ----------
cat("\n########## 【汇总】细胞水平 vs 供体水平 ##########\n")

homc_vs_htc <- if (!is.null(pair_out) && "HomC vs HTC" %in% pair_out$contrast) {
  r <- pair_out[pair_out$contrast == "HomC vs HTC", ]
  sprintf("配对均值差 %.4f, P = %.3g, dz = %.2f (n = %d 供体)",
          r$mean_diff, r$P_paired_t, r$cohens_dz, r$n_donors)
} else "NA"

cmp <- data.frame(
  quantity = c("PIEZO1 亚型差异 (7 亚型 omnibus)",
               "PIEZO1 vs 分化轴 (趋势)",
               "Q1稳态 -> Q4肥大 差异",
               "HomC vs HTC"),
  cell_level = c(
    sprintf("Kruskal-Wallis P = %.3g (n = %d 细胞)", kw_cell$p.value, length(pz_cell)),
    sprintf("rho = %.3f, P = %.3g (n = %d 细胞)", ct_cell$estimate, ct_cell$p.value, length(pz_cell)),
    sprintf("%.3f vs %.3f (差 %.3f)", prof_cell[["Q1_homeo"]], prof_cell[["Q4_hyp"]],
            prof_cell[["Q1_homeo"]] - prof_cell[["Q4_hyp"]]),
    sprintf("均值差 %.3f", mean(pz_cell[md$subtype == "HomC"]) - mean(pz_cell[md$subtype == "HTC"]))),
  donor_level = c(
    if (!is.na(fr_out$P))
      sprintf("Friedman P = %.3g (n = %d 供体, Kendall W = %.2f)",
              fr_out$P, fr_out$n_donors, fr_out$kendall_W) else "NA",
    sprintf("供体 rho 均值 %.3f, 单样本 t P = %.3g, 符号检验 %d/%d P = %.3g",
            mean(donor_tab$rho), t_rho$p.value, n_neg, n_tot, sign_p),
    sprintf("配对 t 均值差 %.4f, P = %.3g, 95%%CI [%.3f, %.3f] (n = %d)",
            t_q$estimate, t_q$p.value, t_q$conf.int[1], t_q$conf.int[2], sum(ok)),
    homc_vs_htc),
  stringsAsFactors = FALSE)
print(cmp, row.names = FALSE)
fwrite(cmp, file.path(RES, "40_P0_cell_vs_donor_comparison.csv"))

## ---------- 7. 供图 ----------
dfp <- donor_tab[!is.na(donor_tab$rho), ]
p1 <- ggplot(dfp, aes(x = reorder(patient, rho), y = rho, fill = rho < 0)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"), guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_hline(yintercept = mean(dfp$rho), linetype = "dashed", colour = "#333333", linewidth = 0.5) +
  labs(x = "Donor", y = "Spearman rho (PIEZO1 vs hypertrophic–homeostatic axis)",
       title = "Donor-level association: the correct inferential unit",
       subtitle = sprintf("GSE104782 | n = %d donors | one-sample t P = %.3g | sign test P = %.3g | mean rho = %.3f",
                          nrow(dfp), t_rho$p.value, sign_p, mean(dfp$rho))) +
  coord_flip()

dq <- donor_tab[ok, ]
p2 <- ggplot(dq, aes(x = mean_Q1_homeo, y = mean_Q4_hyp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_segment(aes(xend = mean_Q4_hyp, yend = mean_Q4_hyp),
               colour = "#B2182B", alpha = 0.6, linewidth = 0.4) +
  geom_point(size = 3, colour = "#4C9BD8") +
  coord_equal() +
  labs(x = "Mean PIEZO1, Q1 (most homeostatic)", y = "Mean PIEZO1, Q4 (most hypertrophic)",
       title = "Within-donor Q1 vs Q4",
       subtitle = sprintf("Paired t: mean diff = %.4f, P = %.3g, n = %d donors",
                          t_q$estimate, t_q$p.value, nrow(dq)))

ggsave(file.path(FIG, "fig60_P0_donorlevel_axis.png"), p1 + p2,
       width = 12.5, height = 4.6, dpi = 300)

## ---------- 8. 结论 ----------
cat("\n================================================================\n")
cat(" 结论摘要（用于改写 Methods / Results / Limitations）\n")
cat("================================================================\n")
cat(sprintf("1. 细胞水平 rho = %.3f (P = %.3g)；供体水平 rho 均值 = %.3f ",
            ct_cell$estimate, ct_cell$p.value, mean(donor_tab$rho)))
cat(sprintf("(单样本 t P = %.3g, 符号检验 P = %.3g)。方向一致，\n", t_rho$p.value, sign_p))
cat(sprintf("   但正确推断单元是 %d 个供体而非 %d 个细胞。\n", n_tot, length(pz_cell)))
cat(sprintf("2. 供体水平 rho 95%%CI [%.3f, %.3f]，对应约 %.1f%% 方差解释率，属弱效应。\n",
            t_rho$conf.int[1], t_rho$conf.int[2], 100 * mean(donor_tab$rho)^2))
cat("3. 细胞水平 P 值应降级为描述性/探索性；供体水平检验进入主文本。\n")
cat("4. 若供体水平检验未达显著，稿件须改为‘方向一致但未达供体级显著’，\n")
cat("   不得以细胞水平 P = 5.3e-8 作为头条证据。\n")
cat("================================================================\n")
cat("完成: results/40_* , figs/fig60_P0_donorlevel_axis.png\n")

sink()
