# ============================================================
# 08b_GSE46750_inflammation_gradient.R
# GSE46750: 同一 OA 患者滑膜 炎症区(I) vs 正常/反应区(N/R) 配对 (12 对供体)
# 注意: 培养的滑膜细胞 -> 免疫 marker 解释为"免疫应答基因程序"而非浸润
# ============================================================
suppressMessages({
  library(GEOquery); library(GSVA); library(ggplot2); library(data.table)
})
source("scripts/_bootstrap.R")
outdir <- "results"; figdir <- "figs"

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

## ---------- 1. 读表达矩阵 ----------
g <- getGEO(filename = "data/GSE46750_series_matrix.txt.gz", getGPL = FALSE)
e <- as.matrix(exprs(g)); pd <- pData(g)
ttl <- as.character(pd$title)
pd$region <- ifelse(grepl("inflammatory", ttl), "Inflam", "NormalReact")
pd$subject <- sub(".*subject:\\s*", "", as.character(pd$`characteristics_ch1`))
cat("region:", paste(names(table(pd$region)), table(pd$region), sep="=", collapse=", "),
    "| 供体:", length(unique(pd$subject)), "\n")

## ---------- 2. 解析 GPL10558 注释 ----------
ln <- readLines(gzfile("data/GPL10558.annot.gz"))
i0 <- grep("platform_table_begin", ln)[1]; i1 <- grep("platform_table_end", ln)[1]
hdr <- strsplit(ln[i0 + 1], "\t")[[1]]
sym_col <- grep("symbol", hdr, ignore.case = TRUE)[1]
id_col  <- which(hdr == "ID")[1]; if (is.na(id_col)) id_col <- 1
cat("注释表行数:", i1 - i0 - 2, "| symbol 列:", hdr[sym_col], "\n")
ann <- fread(text = paste(ln[(i0 + 2):(i1 - 1)], collapse = "\n"),
             sep = "\t", header = FALSE, quote = "", fill = TRUE,
             select = c(id_col, sym_col), col.names = c("probe","symbol"))
ann <- ann[symbol != "" & !is.na(symbol)]

## ---------- 3. 最大均值探针合并 ----------
# 若为原始荧光强度 (非 log2) 则变换
if (max(e, na.rm = TRUE) > 100) {
  e[e < 1] <- 1; e <- log2(e)
  cat("原始强度 -> log2 变换完成\n")
}
map <- ann[probe %in% rownames(e)]
e2 <- e[map$probe, , drop = FALSE]
map$mu <- rowMeans(e2, na.rm = TRUE)
map <- map[order(-map$mu)]
map <- map[!duplicated(map$symbol)]
eg <- e[map$probe, , drop = FALSE]; rownames(eg) <- map$symbol
cat("基因级矩阵:", dim(eg)[1], "基因\n")

## ---------- 4. PIEZO1 炎症区 vs 非炎症区 配对 ----------
if ("PIEZO1" %in% rownames(eg)) {
  pz <- as.numeric(eg["PIEZO1", ])
  wide <- data.frame(subject = pd$subject, region = pd$region, pz = pz)
  w <- reshape(wide, idvar = "subject", timevar = "region", direction = "wide")
  names(w) <- sub("^pz\\.", "", names(w))
  ok <- complete.cases(w[, c("Inflam","NormalReact")])
  tt <- t.test(w$Inflam[ok], w$NormalReact[ok], paired = TRUE)
  cat(sprintf("\n[PIEZO1] 炎症区 vs 非炎症区(配对 n=%d): %+0.3f (95%%CI %+.3f~%+.3f) P=%.3g\n",
              sum(ok), tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value))
  # fig27 配对折线
  pd_plot <- data.frame(subject = w$subject[ok], Inflam = w$Inflam[ok], NR = w$NormalReact[ok])
  pl <- rbind(data.frame(subject = pd_plot$subject, region = "Inflam", pz = pd_plot$Inflam),
              data.frame(subject = pd_plot$subject, region = "NormalReact", pz = pd_plot$NR))
  p1 <- ggplot(pl, aes(region, pz, group = subject)) +
    geom_line(color = "grey60", linewidth = 0.5) +
    geom_point(aes(color = region), size = 2.6) +
    scale_color_manual(values = c("Inflam" = "#E15759", "NormalReact" = "#4C9BD8")) +
    labs(x = NULL, y = "PIEZO1 expression (log2)",
         title = "PIEZO1: inflamed vs non-inflamed synovium in same OA patient (paired)",
         subtitle = sprintf("Paired n=%d, paired t-test P=%.3g", sum(ok), tt$p.value)) +
    theme_bw(base_size = 13) + theme(legend.position = "none",
                                      panel.grid.minor = element_blank())
  ggsave(file.path(figdir, "fig27_GSE46750_piezo1_inflam_pairs.png"), p1,
         width = 6.5, height = 5.5, dpi = 300)
} else {
  pz <- rep(NA_real_, ncol(eg)); cat("\nPIEZO1 不在矩阵中!\n")
}

## ---------- 5. ssGSEA 免疫应答程序 + 配对比较 ----------
keep <- sapply(immune_sets, function(g) sum(g %in% rownames(eg)) >= 3)
gs <- immune_sets[keep]
cat("纳入细胞类:", length(gs), "/", length(immune_sets), "\n")
param <- ssgseaParam(exprData = as.matrix(eg), geneSets = gs, normalize = TRUE)
ss <- t(suppressWarnings(gsva(param)))
fwrite(data.frame(sample = rownames(ss), subject = pd$subject, region = pd$region, ss,
                  check.names = FALSE),
       file.path(outdir, "08b_immune_scores_GSE46750.csv"))

tab <- do.call(rbind, lapply(colnames(ss), function(cc) {
  wide <- data.frame(subject = pd$subject, region = pd$region, v = ss[, cc])
  w <- reshape(wide, idvar = "subject", timevar = "region", direction = "wide")
  names(w) <- sub("^v\\.", "", names(w))
  ok <- complete.cases(w[, c("Inflam","NormalReact")])
  wt <- suppressWarnings(wilcox.test(w$Inflam[ok], w$NormalReact[ok], paired = TRUE))
  ct <- suppressWarnings(cor.test(pz, ss[, cc], method = "spearman"))
  data.frame(cell = cc, diff_I_NR = mean(w$Inflam[ok]) - mean(w$NormalReact[ok]),
             P_paired = wt$p.value, rho_all = unname(ct$estimate), P_corr = ct$p.value)
}))
tab$FDR_paired <- p.adjust(tab$P_paired, "BH")
tab$FDR_corr <- p.adjust(tab$P_corr, "BH")
fwrite(tab, file.path(outdir, "08b_GSE46750_inflam_paired_immune.csv"))
cat("配对 FDR<0.05:", ifelse(any(tab$FDR_paired < 0.05),
    paste(tab$cell[tab$FDR_paired < 0.05], collapse = ", "), "无"),
    "\nPIEZO1 相关 FDR<0.05:", ifelse(any(tab$FDR_corr < 0.05),
    paste(tab$cell[tab$FDR_corr < 0.05], collapse = ", "), "无"), "\n")

## ---------- 6. 图 ----------
theme_set(theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank()))
# fig28 免疫应答 配对箱线
sl <- data.frame(subject = rep(pd$subject, ncol(ss)), region = rep(pd$region, ncol(ss)),
                 cell = rep(colnames(ss), each = nrow(ss)),
                 score = as.numeric(ss))
sigd <- tab$cell[tab$FDR_paired < 0.05]
p2 <- ggplot(sl, aes(region, score, fill = region)) +
  geom_boxplot(alpha = 0.75, width = 0.55, outlier.shape = NA) +
  facet_wrap(~ cell, scales = "free_y") +
  scale_fill_manual(values = c("Inflam" = "#E15759", "NormalReact" = "#4C9BD8")) +
  labs(x = NULL, y = "ssGSEA score",
       title = "Synovial immune-response programs: inflamed vs non-inflamed (same patient, GSE46750)",
       subtitle = paste0("Cultured synovial cells; FDR<0.05: ",
                         ifelse(length(sigd), paste(sigd, collapse=", "), "none"))) +
  theme(legend.position = "none", strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(figdir, "fig28_GSE46750_inflam_immune.png"), p2,
       width = 13, height = 9, dpi = 300)

# fig29 PIEZO1 x 免疫应答 相关散点 (Top4)
topc <- tab$cell[order(-abs(tab$rho_all))][1:4]
pdsc <- data.frame(pz = pz, region = pd$region,
                   stack(as.data.frame(ss[, topc, drop = FALSE])),
                   check.names = FALSE)
names(pdsc)[3:4] <- c("score","cell")
p3 <- ggplot(pdsc, aes(pz, score)) +
  geom_point(aes(color = region), size = 2.4) +
  geom_smooth(method = "lm", se = TRUE, color = "grey35", linewidth = 0.6) +
  facet_wrap(~ cell, scales = "free") +
  scale_color_manual(values = c("Inflam" = "#E15759", "NormalReact" = "#4C9BD8")) +
  labs(x = "PIEZO1 expression (log2)", y = "ssGSEA score",
       title = "PIEZO1 vs immune-response program correlation (GSE46750, n=24)") +
  theme(legend.position = "top", strip.text = element_text(size = 9))
ggsave(file.path(figdir, "fig29_GSE46750_piezo1_immune_scatter.png"), p3,
       width = 10, height = 7.5, dpi = 300)

cat("\n完成: fig27-fig29 + 08b 结果表\n")
