# 10_scrna_piezo1.R — PIEZO1 单细胞定位分析 (GSE104782, OA 软骨)
# 数据: 24153 基因 x 1600 细胞, 10 位 OA 患者软骨, 官方 7 类亚型注释
# 输出: fig36-fig42 + results/10_*.csv

suppressMessages({
  library(data.table); library(Seurat); library(ggplot2)
  library(patchwork); library(readxl); library(pheatmap); library(dplyr)
  library(clusterProfiler); library(org.Hs.eg.db)
})
set.seed(2026)
source("scripts/_bootstrap.R")
figdir <- "figs"; resdir <- "results"
dir.create(figdir, showWarnings = FALSE); dir.create(resdir, showWarnings = FALSE)

theme_set(theme_bw(base_size = 11) +
            theme(panel.grid = element_blank(),
                  plot.title = element_text(face = "bold", size = 11)))

## ---------- 1. 读入数据 ----------
m <- fread("data/scrna/GSE104782_allcells_UMI_count.txt.gz", sep = "\t", header = TRUE)
genes <- m$gene
mat <- as(as.matrix(m[, -1]), "CsparseMatrix")
rownames(mat) <- genes; colnames(mat) <- colnames(m)[-1]
rm(m)

xl <- as.data.frame(read_excel("data/scrna/GSE104782_Table_Cell_quality_information_and_clustering_information.xlsx"))
meta <- xl[match(colnames(mat), xl$Cell), ]
stopifnot(all(meta$Cell == colnames(mat)))
meta$subtype <- meta$Cluster
meta$patient <- sub("_S.*", "", meta$Cell)
meta$s_idx  <- as.integer(sub(".*_S", "", sub("[.].*$", "", meta$Cell)))

cat(sprintf("细胞 %d | 基因 %d | 亚型: %s\n", ncol(mat), nrow(mat),
            paste(names(table(meta$subtype)), collapse = ", ")))

so <- CreateSeuratObject(counts = mat, meta.data = meta[, c("subtype","patient","s_idx",
                                                           "Gene number","Transcript number")])
colnames(so[[]])[4:5] <- c("nGene_paper","nUMI_paper")

## ---------- 2. QC 与过滤 ----------
so[["mt_pct"]] <- PercentageFeatureSet(so, pattern = "^MT-")
keep <- so$subtype[!is.na(so$subtype)]
cells_keep <- colnames(so)[!is.na(so$subtype)]
cat(sprintf("官方注释细胞: %d (剔除未注释 %d)\n", length(cells_keep), ncol(so)-length(cells_keep)))
so <- subset(so, cells = cells_keep)

p_qc <- VlnPlot(so, features = c("nCount_RNA","nFeature_RNA","mt_pct"), ncol = 3,
                group.by = "subtype", pt.size = 0.1) & NoLegend() &
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figdir, "fig36_scrna_QC.png"), p_qc, width = 10, height = 4, dpi = 300)

## ---------- 3. 标准流程 ----------
so <- NormalizeData(so, verbose = FALSE)
so <- FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
so <- ScaleData(so, verbose = FALSE)
so <- RunPCA(so, npcs = 30, verbose = FALSE)
so <- RunUMAP(so, dims = 1:15, verbose = FALSE)
so <- FindNeighbors(so, dims = 1:15, verbose = FALSE)
so <- FindClusters(so, resolution = 0.4, verbose = FALSE)

p1 <- DimPlot(so, group.by = "subtype", label = TRUE, repel = TRUE) +
  ggtitle("OA chondrocyte subtypes (official annotation)")
p2 <- DimPlot(so, group.by = "patient") + ggtitle("Patient origin")
ggsave(file.path(figdir, "fig37_scrna_UMAP_subtypes.png"), p1 + p2,
       width = 12, height = 5, dpi = 300)

## 亚型 marker 确认 (各亚型 Top5)
mk <- FindAllMarkers(so, group.by = "subtype", only.pos = TRUE,
                    logfc.threshold = 0.5, max.cells.per.ident = 200, verbose = FALSE)
topmk <- as.data.frame(mk %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5))
fwrite(topmk, file.path(resdir, "10_scrna_subtype_markers.csv"))
cat("\n=== 各亚型 Top 标志基因 ===\n")
print(topmk[, c("cluster","gene","avg_log2FC","p_val_adj")], row.names = FALSE)

## ---------- 4. PIEZO1 定位: 主结果 ----------
pz <- FetchData(so, "PIEZO1")[, 1]
cat(sprintf("\nPIEZO1: 检出细胞 %d/%d (%.1f%%)\n", sum(pz > 0), length(pz),
            100*mean(pz > 0)))

subtype_order <- c("HomC","ProC","preHTC","HTC","EC","FC","RegC")
sub_lab <- c(HomC = "HomC\nHomeostasis", ProC = "ProC\nProgenitor", preHTC = "preHTC\nPre-hypertrophic",
             HTC = "HTC\nHypertrophic", EC = "EC\nEffector", FC = "FC\nFibrotic", RegC = "RegC\nRegulatory")
so$sub_lab <- factor(unname(sub_lab[as.character(so$subtype)]), levels = unname(sub_lab[subtype_order]))

p_vio <- VlnPlot(so, "PIEZO1", group.by = "sub_lab", pt.size = 0.25,
                 cols = c("#4C9BD8","#59A47F","#E8B04B","#E15759","#B07AA1","#76B7B2","#F28E2B")) +
  stat_summary(fun = mean, geom = "point", size = 1.4) +
  ggtitle("PIEZO1 expression across chondrocyte subtypes", subtitle = "GSE104782, 1464 OA chondrocytes") +
  theme(axis.title.x = element_blank())
p_feat <- FeaturePlot(so, "PIEZO1", reduction = "umap", order = TRUE) +
  ggtitle("PIEZO1 UMAP distribution")
ggsave(file.path(figdir, "fig38_PIEZO1_violin_umap.png"), p_vio + p_feat,
       width = 13, height = 5.5, dpi = 300)

## 亚型水平统计表
dat <- data.frame(subtype = so$subtype, pz = pz,
                  nUMI = so$nCount_RNA, nGene = so$nFeature_RNA)
smry <- do.call(rbind, lapply(split(dat, dat$subtype), function(d) {
  data.frame(subtype = d$subtype[1], n_cells = nrow(d),
             pct_detect = round(100*mean(d$pz > 0), 1),
             mean_expr = round(mean(d$pz[d$pz > 0]), 2),
             mean_all = round(mean(d$pz), 3),
             median_expr = round(median(d$pz[d$pz > 0]), 2))
}))
kw <- kruskal.test(pz ~ subtype, data = dat)
cat("\n=== PIEZO1 亚型分布 (Kruskal-Wallis P =", format.pval(kw$p.value, 3), ") ===\n")
print(smry[order(-smry$mean_expr), ], row.names = FALSE)
fwrite(cbind(smry, KW_P = kw$p.value), file.path(resdir, "10_PIEZO1_subtype_summary.csv"))

## DotPlot: PIEZO1 + 亚型代表 marker
mk_show <- c("PIEZO1","PRG4","CHI3L1","CHI3L2","COL2A1","ACAN","SOX9",
             "COL10A1","RUNX2","IBSP","SPP1","COL1A1","COL1A2","MMP13","CD74","HLA-DRA")
mk_show <- mk_show[mk_show %in% rownames(so)]
p_dot <- DotPlot(so, features = mk_show, group.by = "sub_lab") +
  coord_flip() +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  ggtitle("Subtype marker genes and PIEZO1") +
  theme(axis.text.x = element_text(angle = 0))
ggsave(file.path(figdir, "fig39_PIEZO1_dotplot.png"), p_dot, width = 7, height = 8, dpi = 300)

## ---------- 5. PIEZO1+ vs PIEZO1- 差异分析 ----------
so$pz_pos <- factor(ifelse(pz > 0, "PIEZO1+", "PIEZO1-"), levels = c("PIEZO1-","PIEZO1+"))
cat("\nPIEZO1+ 细胞:", sum(pz > 0), " PIEZO1-:", sum(pz == 0), "\n")
cat("检出与测序深度的关系: Spearman rho =",
    cor.test(as.numeric(pz > 0), so$nCount_RNA, method = "spearman")$estimate,
    "P =", format.pval(cor.test(as.numeric(pz > 0), so$nCount_RNA, method = "spearman")$p.value, 3), "\n")

Idents(so) <- "pz_pos"
de <- FindMarkers(so, ident.1 = "PIEZO1+", ident.2 = "PIEZO1-",
                  test.use = "wilcox", logfc.threshold = 0.25, verbose = FALSE)
de$gene <- rownames(de); de <- de[order(-de$avg_log2FC), ]
fwrite(de, file.path(resdir, "10_PIEZO1pos_vs_neg_DE.csv"))
cat("\n=== PIEZO1+ 特异上调 Top15 ===\n")
print(head(de[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 15), row.names = FALSE)
cat("=== PIEZO1+ 特异下调 Top10 ===\n")
print(tail(de[, c("gene","avg_log2FC","pct.1","pct.2","p_val_adj")], 10), row.names = FALSE)

## 深度匹配敏感性: 仅用 nUMI 中位数以上细胞
so2 <- subset(so, subset = nCount_RNA > median(so$nCount_RNA))
Idents(so2) <- so2$pz_pos
de2 <- FindMarkers(so2, ident.1 = "PIEZO1+", ident.2 = "PIEZO1-", test.use = "wilcox",
                   logfc.threshold = 0.25, verbose = FALSE)
de2$gene <- rownames(de2)
fwrite(de2, file.path(resdir, "10_PIEZO1pos_vs_neg_DE_depthmatched.csv"))

## GO 富集 (PIEZO1+ 上调基因)
up <- de$gene[de$avg_log2FC > 0.25 & de$p_val_adj < 0.05]
cat("\n富集用上调基因:", length(up), "\n")
if (length(up) >= 10) {
  eg <- bitr(up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  ego <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                  pAdjustMethod = "BH", qvalueCutoff = 0.1, readable = TRUE)
  ego_df <- as.data.frame(ego)
  fwrite(ego_df, file.path(resdir, "10_GO_PIEZO1pos_up.csv"))
  if (nrow(ego_df)) {
    e10 <- head(ego_df[order(ego_df$p.adjust), ], 15)
    e10$Description <- substr(e10$Description, 1, 55)
    p_go <- ggplot(e10, aes(Count, reorder(Description, -p.adjust))) +
      geom_col(aes(fill = p.adjust), width = 0.7) +
      scale_fill_gradient(low = "#B2182B", high = "#92C5DE") +
      labs(y = NULL, x = "Gene count",
           title = "GO enrichment of genes up-regulated in PIEZO1+ cells",
           subtitle = sprintf("%d FDR<0.05 genes", length(up))) +
      theme(axis.text.y = element_text(size = 8))
    ggsave(file.path(figdir, "fig40_GO_PIEZO1pos.png"), p_go, width = 9, height = 6, dpi = 300)
  }
}

## ---------- 6. PIEZO1 与关键程序基因的相关性 ----------
prog <- c("COL2A1","ACAN","SOX9","PRG4",            # 稳态/软骨
          "COL10A1","RUNX2","MMP13","IBSP","SPP1","ALPL",  # 肥大
          "COL1A1","COL1A2","COL3A1",               # 纤维化
          "CDKN2A","CDKN1A","GLB1",                  # 衰老
          "IL6","CXCL8","CXCL12","CCL2",            # SASP
          "VEGFA","ANGPTL4",                        # 血管
          "CAV1","CAMK2D","CALM1",                  # 钙信号
          "MMP1","MMP3","ADAMTS5","TNF","IL1B")
prog <- prog[prog %in% rownames(so)]
expr <- as.matrix(GetAssayData(so, layer = "data")[prog, ])
cor_tab <- do.call(rbind, lapply(rownames(expr), function(g) {
  ct <- suppressWarnings(cor.test(pz, expr[g, ], method = "spearman"))
  data.frame(gene = g, rho = round(unname(ct$estimate), 3),
             P = ct$p.value, FDR = p.adjust(ct$p.value, "BH"))
}))
cor_tab <- cor_tab[order(-abs(cor_tab$rho)), ]
fwrite(cor_tab, file.path(resdir, "10_PIEZO1_program_correlation.csv"))
cat("\n=== PIEZO1 与程序基因 Spearman (Top15) ===\n")
print(head(cor_tab, 15), row.names = FALSE)

sig_cor <- cor_tab[cor_tab$FDR < 0.05 & cor_tab$gene != "PIEZO1", ]
if (nrow(sig_cor)) {
  p_cor <- ggplot(sig_cor, aes(rho, reorder(gene, rho), fill = rho > 0)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("#2166AC","#B2182B"), guide = "none") +
    labs(x = "Spearman rho with PIEZO1", y = NULL,
         title = "PIEZO1-correlated program genes (FDR<0.05)",
         subtitle = "Red = positive, Blue = negative | OA chondrocyte single-cell level") +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(figdir, "fig41_PIEZO1_correlation_bars.png"), p_cor,
         width = 8, height = 0.35*nrow(sig_cor) + 2, dpi = 300)
}

## ---------- 7. 患者/样本一致性 ----------
samp <- data.frame(cell = colnames(so), subtype = so$subtype,
                   patient = so$patient, s_idx = so$s_idx,
                   pz = pz, pz_pos = pz > 0)
samp_smry <- aggregate(cbind(pz, pz_pos) ~ patient + s_idx, data = samp, FUN = mean)
names(samp_smry)[3:4] <- c("mean_pz","frac_detect")
samp_smry$frac_detect <- round(100*samp_smry$frac_detect, 1)
fwrite(samp_smry, file.path(resdir, "10_sample_summary.csv"))

## 患者 x 亚型 平均表达热图
pzm <- tapply(pz, list(samp$patient, samp$subtype), mean)
pzm[is.na(pzm)] <- 0
pzm <- pzm[, subtype_order[subtype_order %in% colnames(pzm)]]
p_hm <- pheatmap(pzm, cluster_cols = FALSE, scale = "none",
                 color = colorRampPalette(c("white","#F4A582","#B2182B"))(100),
                 main = "PIEZO1 mean expression: patient x subtype",
                 filename = file.path(figdir, "fig42_PIEZO1_patient_heatmap.png"),
                 width = 6, height = 5.5, dpi = 300, silent = TRUE)

## S 指数趋势 (探索性)
cat("\n=== PIEZO1 按样本序号 S0-S4 (探索性) ===\n")
print(aggregate(pz ~ s_idx, data = samp, FUN = function(x) round(mean(x), 3)))
cat("注: S0-S4 为同一患者软骨的不同取样部位, 原文未明确分区含义\n")

## 亚型组成 x PIEZO1 高低 (每患者)
cat("\n完成: fig36-fig42 + results/10_*.csv\n")
saveRDS(so, "data/scrna/GSE104782_seurat.rds")
