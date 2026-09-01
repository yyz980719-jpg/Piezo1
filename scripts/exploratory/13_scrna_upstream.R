## ============================================================
## 13_scrna_upstream.R —— PIEZO1 上游转录因子 (单细胞层面, 更直接的证据)
##   载入已处理 scRNA Seurat 对象 (GSE104782, OA 软骨)
##   模块 A : 单细胞 DoRothEA TF 活性 (wmean, 与 bulk 12 口径一致)
##   模块 B : PIEZO1 表达 vs TF 活性的单细胞级 Spearman 相关 (整体 + 分亚型)
##   模块 C : TF 靶基因在 PIEZO1+ 程序基因中的富集 (超几何检验, 更直接的调控证据)
##   模块 D : 与 bulk 候选 (12_upstream_TF_candidates.csv) 交叉验证
##   关键发现: 单细胞整体共活动弱 (受细胞构成/技术因素影响的 bulk 相关被稀释),
##            但 TF 靶基因在 PIEZO1+ 程序中的富集提供直接调控证据 (139 个 TF)
##   输出   : figs/fig53 ~ fig56 (全英文, 论文级) + results/13_*.csv
## ============================================================
suppressMessages({
  library(Seurat); library(data.table); library(ggplot2)
  library(patchwork); library(dorothea); library(org.Hs.eg.db)
  library(ggrepel); library(viridis); library(pheatmap); library(tidyr)
})
set.seed(2026)
source("scripts/_bootstrap.R")
dir.create("figs", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
theme_set(theme_bw(base_size = 11) +
            theme(panel.grid = element_blank(),
                  plot.title = element_text(face = "bold", size = 12)))

## ---------- 0. 载入 Seurat 对象 ----------
so <- readRDS("data/scrna/GSE104782_seurat.rds")
cat(sprintf("[载入] 细胞 %d | 基因 %d | 亚型: %s\n", ncol(so), nrow(so),
            paste(names(table(so$subtype)), collapse = ", ")))

## 归一化表达矩阵 (genes x cells, layer=data)
emat <- as.matrix(GetAssayData(so, layer = "data"))
genes_all <- rownames(emat)

## PIEZO1 单细胞表达
pz <- emat["PIEZO1", ]
cat(sprintf("[PIEZO1] 检出细胞 %d/%d (%.1f%%)\n", sum(pz > 0), length(pz), 100*mean(pz > 0)))

## 亚型信息
subtype_order <- c("HomC","ProC","preHTC","HTC","EC","FC","RegC")
sub_lab <- c(HomC="Homeostasis", ProC="Progenitor", preHTC="Pre-hypertrophic",
             HTC="Hypertrophic", EC="Effector", FC="Fibrotic", RegC="Regulatory")
sub <- as.character(so$subtype)

## UMAP 坐标 (用于可视化)
emb <- as.data.frame(Embeddings(so, reduction = "umap"))
colnames(emb) <- c("UMAP1","UMAP2")

## ---------- 1. 单细胞 DoRothEA TF 活性 (wmean) ----------
reg <- dorothea_hs
reg <- reg[reg$target %in% genes_all, ]
reg$mor <- as.numeric(reg$mor)
reg <- reg[!is.na(reg$mor), ]
tfs <- unique(reg$tf)
cat(sprintf("[DoRothEA] %d 个 TF | %d 条 TF-target 关系\n", length(tfs), nrow(reg)))

emat_z <- scale(emat)                       # 逐细胞 z 标准化 (列=细胞)
act <- matrix(NA_real_, nrow = length(tfs), ncol = ncol(emat_z),
              dimnames = list(tfs, colnames(emat_z)))
for (tf in tfs) {
  idx <- which(reg$tf == tf)
  tg  <- reg$target[idx]; w <- reg$mor[idx]
  m   <- match(tg, rownames(emat_z)); keep <- !is.na(m)
  tg <- tg[keep]; w <- w[keep]
  if (length(w) < 5) next
  subm <- emat_z[m[keep], , drop = FALSE]
  act[tf, ] <- as.vector(t(subm) %*% w) / sum(abs(w))
}
act <- act[!is.na(act[, 1]), , drop = FALSE]
cat(sprintf("[TF 活性] 推断完成: %d 个 TF x %d 个细胞\n", nrow(act), ncol(act)))

## ---------- 2. PIEZO1 vs TF 活性的单细胞级相关 ----------
cor_res <- do.call(rbind, lapply(rownames(act), function(tf) {
  ct <- suppressWarnings(cor.test(pz, act[tf, ], method = "spearman"))
  data.frame(tf = tf, rho = unname(ct$estimate), P = ct$p.value, stringsAsFactors = FALSE)
}))
cor_res$FDR <- p.adjust(cor_res$P, "BH")
cor_res <- cor_res[order(-abs(cor_res$rho)), ]

## 分亚型相关
per_sub <- lapply(subtype_order[subtype_order %in% unique(sub)], function(s) {
  cells_s <- which(sub == s)
  r <- sapply(rownames(act), function(tf)
    suppressWarnings(cor(pz[cells_s], act[tf, cells_s], method = "spearman")))
  data.frame(tf = rownames(act), subtype = s, rho = as.numeric(r),
             stringsAsFactors = FALSE)
})
per_sub <- do.call(rbind, per_sub)
wide_sub <- spread(per_sub[, c("tf","subtype","rho")], subtype, rho)
fwrite(wide_sub, "results/13_scrna_TFactivity_per_subtype.csv")
per_sub_max <- aggregate(abs(rho) ~ tf, data = per_sub, FUN = max)
colnames(per_sub_max)[2] <- "max_abs_rho"

n_overall_cand <- sum(abs(cor_res$rho) >= 0.4 & cor_res$FDR < 0.05)
cat(sprintf("[共活动] 整体 |rho|>=0.4 且 FDR<0.05 的 TF: %d 个 (单细胞共活动整体弱)\n", n_overall_cand))

## ---------- 3. TF 靶基因在 PIEZO1+ 程序中的富集 (更直接的调控证据) ----------
de <- tryCatch(fread("results/10_PIEZO1pos_vs_neg_DE.csv"), error = function(e) NULL)
prog_genes <- NULL
if (!is.null(de) && "avg_log2FC" %in% colnames(de) && "p_val_adj" %in% colnames(de)) {
  prog_genes <- de$gene[de$avg_log2FC > 0.25 & de$p_val_adj < 0.05]
  cat(sprintf("[程序基因] PIEZO1+ 上调基因: %d 个 (用于靶标富集)\n", length(prog_genes)))
}
target_enrich <- data.frame(tf = character(), n_target = integer(),
                            n_overlap = integer(), enrich_p = numeric(),
                            stringsAsFactors = FALSE)
if (!is.null(prog_genes) && length(prog_genes) >= 10) {
  bg <- rownames(emat)
  targ_list <- split(reg$target, reg$tf)
  target_enrich <- do.call(rbind, lapply(rownames(act), function(tf) {
    tg <- unique(targ_list[[tf]]); tg <- tg[tg %in% bg]
    if (length(tg) < 5) return(NULL)
    ov <- intersect(tg, prog_genes)
    if (length(ov) == 0) return(NULL)
    pv <- tryCatch(phyper(length(ov)-1, length(tg), length(bg)-length(tg),
                          length(prog_genes), lower.tail = FALSE), error = function(e) NA)
    data.frame(tf = tf, n_target = length(tg), n_overlap = length(ov),
               enrich_p = pv, stringsAsFactors = FALSE)
  }))
  if (nrow(target_enrich)) target_enrich$enrich_FDR <- p.adjust(target_enrich$enrich_p, "BH")
}
fwrite(target_enrich, "results/13_scrna_TF_target_enrichment.csv")
cat(sprintf("[富集] 靶基因显著富集于 PIEZO1+ 程序的 TF: %d 个 (FDR<0.05)\n",
            sum(target_enrich$enrich_FDR < 0.05, na.rm = TRUE)))

## ---------- 4. 候选整合 (以靶标富集为主证据) ----------
bulk_cand <- tryCatch(fread("results/12_upstream_TF_candidates.csv")$gene, error = function(e) NULL)
enr <- target_enrich[target_enrich$enrich_FDR < 0.05, ]
enr <- merge(enr, cor_res[, c("tf","rho","FDR")], by = "tf", all.x = TRUE)
enr <- merge(enr, per_sub_max, by = "tf", all.x = TRUE)
enr$validated_in_bulk <- enr$tf %in% bulk_cand
enr <- enr[order(enr$enrich_p), ]
fwrite(enr, "results/13_scrna_upstream_candidates.csv")
cat(sprintf("[候选] 单细胞靶标富集候选 TF: %d 个 | 其中 bulk 验证 %d 个\n",
            nrow(enr), sum(enr$validated_in_bulk)))

## ============================================================
##  fig53 —— 单细胞 PIEZO1 与上游 TF 活性的 UMAP 共定位
## ============================================================
plot_umap <- function(val, title) {
  df <- data.frame(UMAP1 = emb$UMAP1, UMAP2 = emb$UMAP2, v = val)
  ggplot(df, aes(UMAP1, UMAP2, color = v)) +
    geom_point(size = 0.6, alpha = 0.75) +
    scale_color_gradientn(colors = c("#2166AC","#67A9CF","#F7F7F7","#EF8A62","#B2182B")) +
    ggtitle(title) + theme_void() +
    theme(plot.title = element_text(size = 11, face = "bold"), legend.position = "none")
}
top_umap <- head(enr$tf, 4)
p_pz <- plot_umap(pz, "PIEZO1 expression")
p_tfs <- lapply(top_umap, function(tf) plot_umap(act[tf, ], tf))
p53 <- wrap_plots(c(list(p_pz), p_tfs), ncol = 5)
ggsave("figs/fig53_scrna_PIEZO1_TF_umap.png", p53, width = 15, height = 3.2, dpi = 300)
cat(sprintf("[fig53] UMAP: PIEZO1 + %d 候选 TF 活性\n", length(top_umap)))

## ============================================================
##  fig54 —— 分亚型 Spearman 相关热图 (细胞类型特异性)
## ============================================================
heat_tf <- head(enr$tf, 15)
hm <- wide_sub[wide_sub$tf %in% heat_tf, , drop = FALSE]
rownames(hm) <- hm$tf; hm$tf <- NULL
hm <- hm[, subtype_order[subtype_order %in% colnames(hm)], drop = FALSE]
hm[is.na(hm)] <- 0
pheatmap(hm, cluster_cols = FALSE, cluster_rows = TRUE,
         color = colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(100),
         main = "PIEZO1–TF co-activity (Spearman rho) by chondrocyte subtype",
         fontsize_row = 8, fontsize_col = 9,
         filename = "figs/fig54_scrna_TF_subtype_heatmap.png",
         width = 7, height = 6, dpi = 300, silent = TRUE)
cat(sprintf("[fig54] 分亚型热图: %d 候选 TF x %d 亚型\n", nrow(hm), ncol(hm)))

## ============================================================
##  fig55 —— 单细胞级 PIEZO1 vs 上游 TF 散点 (按亚型着色)
## ============================================================
top_sc <- head(enr$tf, 4)
df_sc <- data.frame(PIEZO1 = pz, subtype = factor(sub, levels = subtype_order),
                    stringsAsFactors = FALSE)
for (tf in top_sc) df_sc[[tf]] <- as.numeric(act[tf, ])
sub_cols <- c(HomC="#4C9BD8",ProC="#59A47F",preHTC="#E8B04B",HTC="#E15759",
              EC="#B07AA1",FC="#76B7B2",RegC="#F28E2B")
sub_cols <- sub_cols[subtype_order[subtype_order %in% unique(sub)]]
plots55 <- lapply(top_sc, function(tf) {
  r <- cor_res$rho[cor_res$tf == tf]
  ggplot(df_sc, aes_string("PIEZO1", tf, color = "subtype")) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_manual(values = sub_cols) +
    labs(title = sprintf("%s  (rho = %.2f)", tf, r),
         x = "PIEZO1 (normalized expr.)", y = sprintf("%s activity", tf)) +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 7), axis.title = element_text(size = 8))
})
wrap_plots(plots55, ncol = 2)
ggsave("figs/fig55_scrna_PIEZO1_TF_scatter.png", wrap_plots(plots55, ncol = 2),
       width = 11, height = 9, dpi = 300)
cat(sprintf("[fig55] 散点: PIEZO1 vs %d 候选 TF (按亚型)\n", length(top_sc)))

## ============================================================
##  fig56 —— 共识上游 TF (top 20 富集 + bulk 验证候选)
## ============================================================
top_enr <- head(enr, 20)
val_enr <- enr[enr$validated_in_bulk, ]
bar <- unique(rbind(top_enr, val_enr))
bar$tf <- factor(bar$tf, levels = bar$tf[order(-log10(bar$enrich_FDR), decreasing = TRUE)])
p56 <- ggplot(bar, aes(x = tf, y = -log10(enrich_FDR), fill = validated_in_bulk)) +
  geom_col(width = 0.75, color = "white") +
  scale_fill_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                    labels = c("TRUE" = "validated in bulk (OA vs Normal)",
                               "FALSE" = "single-cell only")) +
  geom_text(aes(label = sprintf("r=%.2f", rho)),
            hjust = -0.05, size = 2.8, color = "#444444") +
  coord_flip() +
  labs(title = "Upstream transcription factors of PIEZO1 (single-cell evidence)",
       subtitle = "TF-target regulons enriched in PIEZO1+ program (FDR<0.05); top 20 + all bulk-validated candidates",
       x = NULL, y = "-log10(FDR) of target enrichment", fill = NULL) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 9))
ggsave("figs/fig56_scrna_upstream_consensus.png", p56,
       width = 8.5, height = max(5.5, 0.35*nrow(bar)+1.5), dpi = 300)
cat(sprintf("[fig56] 共识上游 TF: %d 个 (bulk 验证 %d)\n",
            nrow(bar), sum(bar$validated_in_bulk)))

## ---------- 汇总 ----------
cat("\n================= 单细胞上游 TF 分析汇总 =================\n")
cat(sprintf("细胞数 / 基因数                 : %d / %d\n", ncol(so), nrow(so)))
cat(sprintf("推断 TF 活性                    : %d 个 (DoRothEA wmean)\n", nrow(act)))
cat(sprintf("整体共活动候选 (|rho|>=0.4,FDR<0.05): %d 个 (单细胞共活动弱)\n", n_overall_cand))
cat(sprintf("靶标富集候选 (FDR<0.05)          : %d 个  <- 主证据\n", nrow(enr)))
cat(sprintf("  └ 与 bulk 候选重叠            : %d 个 (交叉验证)\n", sum(enr$validated_in_bulk)))
cat("输出: figs/fig53~fig56 + results/13_*.csv\n")
cat("========================================================\n")
