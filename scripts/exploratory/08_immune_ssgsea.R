# ============================================================
# 08_immune_ssgsea.R  免疫浸润 x PIEZO1 (骨关节炎)
# 数据集:
#   GSE55235  滑膜   10 Healthy vs 10 OA      (GPL96)
#   GSE82107  滑膜    7 Healthy vs 10 OA      (GPL570)
#   GSE51588  软骨下骨 20 OA配对(MT负重/LT非负重) + 5 Normal配对 (GPL13497)
#   GSE12021  滑膜   (若下载成功则纳入, GPL96)  [外部验证]
# 分析:
#   A. ssGSEA 定量免疫细胞浸润(18类经典marker)
#   B. OA vs 对照 组间免疫差异 (Wilcoxon + BH)
#   C. PIEZO1 x 免疫细胞 Spearman 相关 (跨数据集一致性)
#   D. GSE51588 负重/非负重配对免疫比较 + PIEZO1-delta x 免疫-delta 相关
# ============================================================
suppressMessages({
  library(GEOquery); library(GSVA); library(ggplot2); library(data.table)
  library(hgu133a.db); library(hgu133plus2.db); library(pROC)
})
source("scripts/_bootstrap.R")
outdir <- "results"; figdir <- "figs"
dir.create(outdir, showWarnings = FALSE); dir.create(figdir, showWarnings = FALSE)

## ---------- 免疫细胞 marker 集合 (18 类, 经典文献 marker 汇总) ----------
immune_sets <- list(
  "T_cells"          = c("CD3D","CD3E","CD3G","CD2","CD28","TRAC"),
  "CD8_T_cells"      = c("CD8A","CD8B","GZMK","GZMA","GZMB"),
  "CD4_T_cells"      = c("CD4","IL7R","CD3D","LEF1"),
  "Tregs"            = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  "Th1_cells"        = c("IFNG","CXCR3","STAT1","TBX21"),
  "Th2_cells"        = c("GATA3","IL4","IL5","IL13","CCR4"),
  "Th17_cells"       = c("RORC","IL17A","IL23R","CCR6"),
  "B_cells"          = c("CD19","MS4A1","CD79A","CD79B","CD22","TCL1A"),
  "Plasma_cells"     = c("MZB1","DERL3","TNFRSF17","SDC1","XBP2"),
  "NK_cells"         = c("NKG7","KLRD1","GNLY","NCR1","KIR2DL3"),
  "Monocytes"        = c("CD14","LST1","FCN1","S100A12","CD33"),
  "Macrophages_M1"   = c("CD86","SOCS3","IL1B","CXCL9","CXCL10"),
  "Macrophages_M2"   = c("CD163","MRC1","MSR1","CCL18","CD14"),
  "Dendritic_cells"  = c("FCER1A","CLEC10A","CD1C","BATF3"),
  "Plasmacytoid_DC"  = c("CLEC4C","LILRA4","IL3RA","TCF4"),
  "Neutrophils"      = c("FCGR3B","CSF3R","S100A8","S100A9","CXCR1"),
  "Mast_cells"       = c("TPSAB1","TPSB2","CPA3","KIT","MS4A2"),
  "Eosinophils"      = c("CCR3","SIGLEC8","PRG2","EPX")
)

## ---------- 工具函数 ----------
# 最大均值探针合并: 每 symbol 保留平均表达最高的探针
collapse_probes <- function(e, map) {
  map <- map[map$probe %in% rownames(e) & !is.na(map$symbol) & map$symbol != "", ]
  e2  <- e[map$probe, , drop = FALSE]
  mu  <- rowMeans(e2, na.rm = TRUE)
  map$mu <- mu
  map <- map[order(-map$mu), ]
  map <- map[!duplicated(map$symbol), ]
  out <- e[map$probe, , drop = FALSE]
  rownames(out) <- map$symbol
  out
}

# 鲁棒分组: 优先扫描 characteristics_ch* 列 (疾病状态通常在此处), 避免 title/source_name 误判
robust_group <- function(pd) {
  cols <- grep("characteristics", colnames(pd), value = TRUE)
  if (!length(cols)) cols <- colnames(pd)
  for (cc in cols) {
    v <- as.character(pd[[cc]])
    if (any(grepl("osteoarthrit", v, ignore.case = TRUE)) ||
        any(toupper(trimws(v)) == "OA")) {
      vv <- toupper(trimws(sub(".*:", "", v)))
      g <- ifelse(grepl("OSTEOARTH", vv) | vv == "OA", "OA", NA)
      g <- ifelse(is.na(g) & (grepl("HEALTHY", vv) | grepl("NORMAL", vv)), "Healthy", g)
      if (sum(!is.na(g)) >= 5) return(g)
    }
  }
  rep(NA, nrow(pd))
}

# ssGSEA (GSVA 2.x API, 兼容旧版)
run_ssgsea <- function(mat) {
  keep <- sapply(immune_sets, function(g) sum(g %in% rownames(mat)) >= 3)
  gs   <- immune_sets[keep]
  dropped <- names(immune_sets)[!keep]
  if (length(dropped)) cat("  丢弃基因<3的细胞类:", paste(dropped, collapse=", "), "\n")
  param <- ssgseaParam(exprData = as.matrix(mat), geneSets = gs, normalize = TRUE)
  suppressWarnings(gsva(param))
}

# ---------- 读单个 GEO 数据集并得到 gene-level 矩阵 ----------
load_geo <- function(gsm_file, platform) {
  g <- getGEO(filename = gsm_file, getGPL = FALSE)
  e <- exprs(g); pd <- pData(g)
  e <- e[rownames(e) != "", , drop = FALSE]
  if (max(e, na.rm = TRUE) > 100) {
    e[e < 1] <- 1; e <- log2(e)
    cat("原始强度 -> log2 变换 (max > 100)\n")
  }
  if (platform == "GPL96") {
    sy <- toTable(hgu133aSYMBOL); map <- data.frame(probe = sy$probe_id, symbol = sy$symbol)
  } else if (platform == "GPL570") {
    sy <- toTable(hgu133plus2SYMBOL); map <- data.frame(probe = sy$probe_id, symbol = sy$symbol)
  } else {  # GPL13497
    map <- fread("data/GPL13497_annot_full.tsv", header = FALSE,
                 col.names = c("probe","symbol"), na.strings = c("", "NA", "---"))
    map <- unique(as.data.frame(map))
  }
  eg <- collapse_probes(e, map)
  list(e = eg, pd = pd)
}

## ============================================================
## 数据集 1/2/4: 滑膜 (OA vs Healthy)
## ============================================================
datasets <- list()

# --- GSE55235 (GPL96) ---
cat("\n===== GSE55235 滑膜 =====\n")
r <- load_geo("data/GSE55235_series_matrix.txt.gz", "GPL96")
grp <- robust_group(r$pd); r$pd$group <- grp
cat("分组:", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "), "\n")
datasets[["GSE55235_Synovium"]] <- list(e = r$e, pd = r$pd)

# --- GSE82107 (GPL570) ---
cat("\n===== GSE82107 滑膜 =====\n")
r <- load_geo("data/GSE82107_series_matrix.txt.gz", "GPL570")
grp <- robust_group(r$pd); r$pd$group <- grp
cat("分组:", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "), "\n")
datasets[["GSE82107_Synovium"]] <- list(e = r$e, pd = r$pd)

# --- GSE12021 (GPL96, 外部验证, 若下载成功) ---
f12021 <- "data/GSE12021_series_matrix.txt.gz"
if (file.exists(f12021) && file.size(f12021) > 1e6) {
  cat("\n===== GSE12021 滑膜 (外部验证) =====\n")
  ok <- tryCatch({ gz <- gzfile(f12021); lns <- readLines(gz, n = 5); close(gz); TRUE },
                 error = function(e) FALSE)
  if (ok) {
    r <- load_geo(f12021, "GPL96")
    grp <- robust_group(r$pd); r$pd$group <- grp
    cat("分组:", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "), "\n")
    # 只保留 OA 与 Healthy
    keep <- r$pd$group %in% c("OA","Healthy")
    datasets[["GSE12021_Synovium"]] <- list(e = r$e[, keep, drop = FALSE], pd = r$pd[keep, ])
    cat("纳入 OA+Healthy 样本:", sum(keep), "\n")
  }
}

## ============================================================
## A+B+C. 滑膜数据集: ssGSEA + 组间比较 + PIEZO1 相关
## ============================================================
corr_all <- list(); diff_all <- list(); score_all <- list()

for (ds in names(datasets)) {
  cat("\n--- ssGSEA:", ds, "---\n")
  d <- datasets[[ds]]
  ss <- run_ssgsea(d$e)
  ss <- t(ss)  # samples x cells
  score_all[[ds]] <- data.frame(sample = rownames(ss), group = d$pd$group, ss,
                                check.names = FALSE)
  fwrite(score_all[[ds]], file.path(outdir, paste0("08_immune_scores_", sub("_.*","",ds), ".csv")))

  # B. 组间差异
  grp <- d$pd$group
  cells <- colnames(ss)
  tab <- do.call(rbind, lapply(cells, function(cc) {
    x <- ss[, cc]
    a <- x[grp == "OA"]; b <- x[grp == "Healthy"]
    w <- suppressWarnings(wilcox.test(a, b))
    data.frame(cell = cc, mean_OA = mean(a), mean_Healthy = mean(b),
               diff = mean(a) - mean(b), P = w$p.value)
  }))
  tab$FDR <- p.adjust(tab$P, "BH")
  tab$dataset <- ds
  diff_all[[ds]] <- tab
  sig <- tab[tab$FDR < 0.05, ]
  cat("FDR<0.05 免疫细胞:", ifelse(nrow(sig), paste(sig$cell, collapse=", "), "无"), "\n")

  # C. PIEZO1 相关
  if ("PIEZO1" %in% rownames(d$e)) {
    pz <- as.numeric(d$e["PIEZO1", ])
    tab2 <- do.call(rbind, lapply(cells, function(cc) {
      ct <- suppressWarnings(cor.test(pz, ss[, cc], method = "spearman"))
      data.frame(cell = cc, rho = unname(ct$estimate), P = ct$p.value)
    }))
    tab2$FDR <- p.adjust(tab2$P, "BH")
    tab2$dataset <- ds
    corr_all[[ds]] <- tab2
    sig2 <- tab2[tab2$FDR < 0.05, ]
    cat("PIEZO1 相关 FDR<0.05:", ifelse(nrow(sig2),
        paste0(sig2$cell, "(rho=", sprintf("%.2f", sig2$rho), ")", collapse = ", "), "无"), "\n")
  }
}

diff_tab <- rbindlist(diff_all)
fwrite(diff_tab, file.path(outdir, "08_immune_groupdiff_OA_vs_Healthy.csv"))
corr_tab <- rbindlist(corr_all)
fwrite(corr_tab, file.path(outdir, "08_PIEZO1_immune_correlation.csv"))

## ============================================================
## D. GSE51588 软骨下骨: 配对区域(负重/非负重)免疫浸润
## ============================================================
cat("\n===== GSE51588 软骨下骨 区域免疫 =====\n")
r <- load_geo("data/GSE51588_series_matrix.txt.gz", "GPL13497")
pd <- r$pd; ttl <- as.character(pd$title)
parts <- strsplit(ttl, "-")
pd$dis <- sapply(parts, `[`, 1)
pd$region <- sapply(parts, `[`, 2)
pd$donor <- paste(pd$dis, sapply(parts, `[`, 3), sep = "-")
cat("region:", paste(names(table(pd$region)), table(pd$region), sep = "=", collapse = ", "), "\n")

ss <- run_ssgsea(r$e); ss <- t(ss)
pz <- as.numeric(r$e["PIEZO1", ])

# OA 患者 MT vs LT 配对比较
oa <- pd$dis == "OA" & pd$region %in% c("MT","LT")
ss_oa <- ss[oa, ]; pd_oa <- pd[oa, ]
regions <- data.frame(sample = rownames(ss_oa), region = pd_oa$region,
                      donor = pd_oa$donor, pz = pz[oa], ss_oa, check.names = FALSE)
fwrite(regions, file.path(outdir, "08_immune_scores_GSE51588_regions.csv"))

tabD <- do.call(rbind, lapply(colnames(ss), function(cc) {
  # 供体配对 Wilcoxon
  wide <- reshape(regions[, c("donor","region",cc)], idvar = "donor",
                 timevar = "region", direction = "wide")
  names(wide) <- sub("^.*\\.", "", names(wide))
  ok <- complete.cases(wide[, c("MT","LT")])
  w <- suppressWarnings(wilcox.test(wide$MT[ok], wide$LT[ok], paired = TRUE))
  # PIEZO1 delta x 免疫 delta 相关
  dw <- data.frame(donor = wide$donor[ok], dlt = wide$LT[ok] - wide$MT[ok])
  pzv <- tapply(pz[oa][pd_oa$region == "LT"], pd_oa$donor[pd_oa$region == "LT"], identity)
  dlt_pz <- pz[oa][pd_oa$region == "LT"] - pz[oa][pd_oa$region == "MT"]
  names(dlt_pz) <- pd_oa$donor[pd_oa$region == "LT"]
  dlt_pz <- dlt_pz[dw$donor]; dw$dz <- dlt_pz
  ct <- suppressWarnings(cor.test(dw$dz, dw$dlt, method = "spearman"))
  data.frame(cell = cc, mean_MT = mean(wide$MT[ok]), mean_LT = mean(wide$LT[ok]),
             diff_LT_MT = mean(wide$LT[ok]) - mean(wide$MT[ok]),
             P_paired = w$p.value, rho_delta = unname(ct$estimate), P_delta = ct$p.value)
}))
tabD$FDR_paired <- p.adjust(tabD$P_paired, "BH")
tabD$FDR_delta  <- p.adjust(tabD$P_delta, "BH")
fwrite(tabD, file.path(outdir, "08_GSE51588_region_paired_immune.csv"))
cat("配对 FDR<0.05:", ifelse(any(tabD$FDR_paired<0.05),
    paste(tabD$cell[tabD$FDR_paired<0.05], collapse=", "), "无"),
    "| delta相关 FDR<0.05:", ifelse(any(tabD$FDR_delta<0.05),
    paste(tabD$cell[tabD$FDR_delta<0.05], collapse=", "), "无"), "\n")

## ============================================================
## 图
## ============================================================
theme_set(theme_bw(base_size = 12) +
          theme(panel.grid.minor = element_blank()))

# fig22: 各滑膜数据集 OA vs Healthy 免疫细胞差异热图 (mean_OA - mean_Healthy)
sc2 <- data.table(diff_tab[, c("dataset","cell","diff","P")])
sc2[, FDR := p.adjust(P, "BH"), by = dataset]
dmat <- dcast(sc2, dataset ~ cell, value.var = "diff")
mat22 <- as.matrix(dmat[, -1]); rownames(mat22) <- dmat$dataset
lab22 <- matrix("", nrow = nrow(mat22), ncol = ncol(mat22))
for (i in seq_len(nrow(sc2))) {
  r <- which(rownames(mat22) == sc2$dataset[i])
  c <- which(colnames(mat22) == sc2$cell[i])
  if (sc2$FDR[i] < 0.05) lab22[r, c] <- "*"
}
png(file.path(figdir, "fig22_immune_diff_heatmap.png"), width = 2200, height = 700, res = 200)
pheatmap::pheatmap(mat22, cluster_rows = FALSE, cluster_cols = TRUE,
                   display_numbers = lab22, number_color = "black",
                   breaks = seq(-0.45, 0.45, length.out = 100),
                   color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
                   main = "OA vs Healthy immune cell infiltration difference (ssGSEA mean difference; * FDR<0.05)",
                   fontsize_row = 11, fontsize_col = 9, angle_col = 45)
dev.off()

# fig23: PIEZO1 x 免疫细胞 相关性热图 (dataset x cell)
rho_mat <- matrix(NA, nrow = length(corr_all), ncol = length(immune_sets),
                 dimnames = list(names(corr_all), names(immune_sets)))
for (ds in names(corr_all)) {
  t <- corr_tab[corr_tab$dataset == ds, ]
  rho_mat[ds, t$cell] <- t$rho
}
pmap <- which(!is.na(rho_mat), arr.ind = TRUE)
lab <- matrix("", nrow = nrow(rho_mat), ncol = ncol(rho_mat))
star <- corr_tab$FDR < 0.05
for (k in seq_len(nrow(pmap))) {
  i <- pmap[k,1]; j <- pmap[k,2]
  t <- corr_tab[corr_tab$dataset == rownames(rho_mat)[i] & corr_tab$cell == colnames(rho_mat)[j], ]
  if (nrow(t) && t$FDR < 0.05) lab[i, j] <- "*"
}
png(file.path(figdir, "fig23_piezo1_immune_corr_heatmap.png"), width = 2400, height = 900, res = 220)
pheatmap::pheatmap(rho_mat, cluster_rows = FALSE, cluster_cols = FALSE,
                   display_numbers = lab, number_color = "black",
                   fontsize_number = 14, breaks = seq(-0.9, 0.9, length.out = 100),
                   color = colorRampPalette(c("#2166AC","white","#B2182B"))(100),
                   main = "PIEZO1 x immune cell infiltration Spearman rho (* FDR<0.05)",
                   fontsize_row = 11, fontsize_col = 10, angle_col = 45)
dev.off()

# fig24: 跨数据集一致性 Top 细胞散点 (在>=2个数据集 |rho|>0.4 且方向一致)
cons <- table(corr_tab$cell[abs(corr_tab$rho) > 0.4])
top_cells <- names(cons[cons >= 2])
if (!length(top_cells)) top_cells <- corr_tab$cell[order(-abs(corr_tab$rho))][1:4]
cat("一致性 Top 细胞:", paste(top_cells, collapse = ", "), "\n")

# fig24: 每个数据集 PIEZO1 vs top 细胞散点
pdat <- rbindlist(lapply(names(score_all), function(ds) {
  x <- score_all[[ds]]; d <- datasets[[ds]]
  if (!"PIEZO1" %in% rownames(d$e)) return(NULL)
  data.frame(dataset = ds, pz = as.numeric(d$e["PIEZO1", ]),
             stack(x[, top_cells, drop = FALSE]), check.names = FALSE)
}))
names(pdat)[3:4] <- c("score","cell")
if (!is.null(pdat) && nrow(pdat)) {
  p <- ggplot(pdat, aes(pz, score)) +
    geom_point(aes(color = dataset), size = 2) +
    geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6) +
    facet_wrap(~ cell, scales = "free") +
    labs(x = "PIEZO1 expression (log2)", y = "ssGSEA score",
         title = "PIEZO1 vs immune cell infiltration correlation (synovium)") +
    theme(legend.position = "top")
  ggsave(file.path(figdir, "fig24_piezo1_immune_scatter.png"), p,
         width = 11, height = 8, dpi = 300)
}

# fig25: GSE51588 区域配对免疫比较 (OA 内 MT vs LT)
reg_long <- data.frame(donor = regions$donor, region = regions$region,
                       stack(regions[, colnames(ss)]), check.names = FALSE)
names(reg_long)[3:4] <- c("score","cell")
sigd <- tabD$cell[tabD$FDR_paired < 0.05]
p <- ggplot(reg_long, aes(region, score, fill = region)) +
  geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA) +
  geom_line(aes(group = donor), color = "grey55", linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(group = donor), size = 1.3, position = position_dodge(0.05)) +
  facet_wrap(~ cell, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("MT" = "#E15759", "LT" = "#4C9BD8")) +
  labs(x = NULL, y = "ssGSEA score",
       title = "Subchondral bone immune infiltration: weight-bearing (MT) vs non-weight-bearing (LT) - OA patient paired",
       subtitle = paste0("Red lines connect same donor; FDR<0.05: ",
                         ifelse(length(sigd), paste(sigd, collapse=", "), "none"))) +
  theme(legend.position = "none", strip.text = element_text(size = 8))
ggsave(file.path(figdir, "fig25_GSE51588_region_immune.png"), p,
       width = 14, height = 10, dpi = 300)

# fig26: PIEZO1-delta x 免疫-delta 散点 (Top 4 by |rho|)
topd <- tabD$cell[order(-abs(tabD$rho_delta))][1:4]
dw_all <- regions[, c("donor","region")]
dw_list <- lapply(topd, function(cc) {
  wide <- reshape(regions[, c("donor","region", cc)], idvar = "donor",
                 timevar = "region", direction = "wide")
  names(wide) <- sub("^.*\\.", "", names(wide))
  ok <- complete.cases(wide[, c("MT","LT")])
  dlt_pz <- pz[oa][pd_oa$region == "LT"] - pz[oa][pd_oa$region == "MT"]
  names(dlt_pz) <- pd_oa$donor[pd_oa$region == "LT"]
  data.frame(dz = dlt_pz[wide$donor[ok]], di = wide$LT[ok] - wide$MT[ok], cell = cc)
})
dw_df <- do.call(rbind, dw_list)
p <- ggplot(dw_df, aes(dz, di)) +
  geom_point(size = 2.2, color = "#845C97") +
  geom_smooth(method = "lm", se = TRUE, color = "#B07AA1", linewidth = 0.7) +
  facet_wrap(~ cell, scales = "free_y") +
  labs(x = "PIEZO1 expression difference (LT - MT, same patient)",
       y = "Immune infiltration difference (LT - MT)",
       title = "Mechanically-related changes: PIEZO1 co-varies with immune cells (subchondral bone)") +
  theme(strip.text = element_text(size = 9))
ggsave(file.path(figdir, "fig26_GSE51588_delta_cor.png"), p,
       width = 9, height = 7, dpi = 300)

cat("\n===== 全部完成 =====\n")
cat("结果表:", list.files(outdir, pattern = "^08_"), "\n")
cat("图:", list.files(figdir, pattern = "^fig2[2-6]"), "\n")
