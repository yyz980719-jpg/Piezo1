## ============================================================
## 12_piezo1_upstream.R —— PIEZO1 上游转录因子 (TF) 分析
##   模块 A : TF 活性推断 (DoRothEA / decoupleR) —— OA vs Normal 差异 TF 活性
##   模块 B : PIEZO1 共表达网络 —— 与 PIEZO1 表达显著相关的候选上游 TF
##   模块 C : 上游候选 TF 富集 + PIEZO1 共表达模块的 GO 功能注释
##   输出   : figs/fig49 ~ fig52 (全英文, 论文级) + results/12_*.csv
## ============================================================
suppressMessages({
  library(GEOquery); library(limma); library(data.table)
  library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db)
  library(igraph); library(viridis); library(pheatmap); library(ggrepel)
})
source("scripts/_bootstrap.R")
dir.create("figs", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

## ---------- 0. 载入 GSE51588 基因级表达矩阵 (与 07 一致) ----------
g   <- getGEO(filename = "data/GSE51588_series_matrix.txt.gz", getGPL = FALSE)
e   <- as.matrix(exprs(g)); pd <- pData(g)
ttl <- as.character(pd$title)
dis <- ifelse(grepl("^OA-", ttl), "OA", "Normal")
reg <- ifelse(grepl("-MT-", ttl), "Medial", "Lateral")
don <- sub("^[^-]+-[^-]+-", "", ttl)
uid <- paste(dis, don, sep = "_")
if (max(e, na.rm = TRUE) > 100) e <- log2(e + 1)
ann <- fread("data/GPL13497_annot_full.tsv", header = FALSE, skip = 1,
             sep = "\t", colClasses = "character",
             col.names = c("probe", "symbol"))
ann <- ann[!is.na(symbol) & nzchar(symbol) &
           grepl("^[A-Za-z][A-Za-z0-9@\\.\\-]*$", symbol)]
common <- intersect(rownames(e), ann$probe)
eg <- e[common, , drop = FALSE]
gl <- ann$symbol[match(common, ann$probe)]
me <- rowMeans(eg, na.rm = TRUE); o <- order(-me)
eg <- eg[o, , drop = FALSE]; gl <- gl[o]
kd <- !duplicated(gl)
emat <- eg[kd, , drop = FALSE]; rownames(emat) <- gl[kd]
storage.mode(emat) <- "numeric"
cat(sprintf("[载入] 样本 %d (OA %d / Normal %d) | 基因 %d | PIEZO1: %s\n",
            ncol(emat), sum(dis=="OA"), sum(dis=="Normal"),
            nrow(emat), "PIEZO1" %in% rownames(emat)))

## ---------- 1. TF 全集 (GO:0003700 DNA 结合转录因子活性) ----------
tf_entrez <- unique(unlist(as.list(org.Hs.egGO2ALLEGS[["GO:0003700"]])))
tf_sym    <- unique(unname(mapIds(org.Hs.eg.db, keys = tf_entrez,
                                  column = "SYMBOL", keytype = "ENTREZID")))
tf_sym    <- tf_sym[!is.na(tf_sym) & tf_sym %in% rownames(emat)]
cat(sprintf("[TF 全集] %d 个转录因子可用于分析\n", length(tf_sym)))

## ---------- 2. PIEZO1 全样本共表达相关性 ----------
pz   <- emat["PIEZO1", ]
rho  <- apply(emat, 1, function(x) suppressWarnings(cor(x, pz, method = "spearman")))
cor_all <- data.frame(gene = names(rho), rho = as.numeric(rho), stringsAsFactors = FALSE)
cor_all <- cor_all[order(-abs(cor_all$rho)), ]
cor_all$is_TF <- cor_all$gene %in% tf_sym
fwrite(cor_all, "results/12_PIEZO1_coexpression.csv")
cat(sprintf("[共表达] 与 PIEZO1 正相关 (rho>=0.5): %d | 负相关 (rho<=-0.5): %d\n",
            sum(cor_all$rho >= 0.5), sum(cor_all$rho <= -0.5)))

## ---------- 3. 模块 A: TF 活性推断 (DoRothEA wmean) ----------
## 实现 DoRothEA 的加权均值 (wmean) 算法: 每个 TF 的活性 =
## 其靶基因表达(按样本 z 标准化)按 TF-靶 权重 (mor) 的加权平均。
use_dorothea <- FALSE
if (requireNamespace("dorothea", quietly = TRUE)) {
  suppressMessages(library(dorothea))
  reg <- dorothea_hs
  reg <- reg[reg$target %in% rownames(emat), ]
  reg$mor <- as.numeric(reg$mor)
  reg <- reg[!is.na(reg$mor), ]
  emat_z <- scale(emat)                      # 逐样本 z 标准化 (列=样本)
  tfs <- unique(reg$tf)
  act <- matrix(NA_real_, nrow = length(tfs), ncol = ncol(emat_z),
                dimnames = list(tfs, colnames(emat_z)))
  for (tf in tfs) {
    idx <- which(reg$tf == tf)
    tg  <- reg$target[idx]; w <- reg$mor[idx]
    m   <- match(tg, rownames(emat_z)); keep <- !is.na(m)
    tg <- tg[keep]; w <- w[keep]
    if (length(w) < 5) next
    sub <- emat_z[m[keep], , drop = FALSE]
    act[tf, ] <- as.vector(t(sub) %*% w) / sum(abs(w))
  }
  act <- act[!is.na(act[, 1]), , drop = FALSE]
  dsg <- model.matrix(~ factor(dis, levels = c("Normal", "OA")))
  fit <- eBayes(lmFit(act, dsg))
  ttA <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  ttA$tf <- rownames(ttA)
  fwrite(ttA, "results/12_TFactivity_OAvsNormal.csv")
  use_dorothea <- TRUE
  cat(sprintf("[DoRothEA wmean] 推断 %d 个 TF 活性 | OA vs Normal 显著差异 %d 个\n",
              nrow(act), sum(ttA$adj.P.Val < 0.05, na.rm = TRUE)))
}

## 降级方案: 若 DoRothEA 不可用, 用 TF 差异表达替代 "TF 活性"
if (!use_dorothea) {
  cat("[降级] 使用 TF 差异表达作为 TF 活性代理\n")
  emTF <- emat[tf_sym, , drop = FALSE]
  dsg  <- model.matrix(~ factor(dis, levels = c("Normal", "OA")))
  fit  <- eBayes(lmFit(emTF, dsg))
  ttA  <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  ttA$tf <- rownames(ttA)
  fwrite(ttA, "results/12_TFactivity_OAvsNormal.csv")
  act <- emTF
}

## ---------- 4. 上游候选 TF 整合 ----------
## 候选 = 与 PIEZO1 强共表达 (|rho|>=0.4) 且 OA 中 TF 活性显著改变 (FDR<0.05)
tfcor_expr <- cor_all[cor_all$is_TF & cor_all$gene %in% ttA$tf, ]
tfcor_expr <- merge(tfcor_expr, ttA[, c("tf","logFC","P.Value","adj.P.Val")],
                    by.x = "gene", by.y = "tf", all.x = TRUE)
tfcor_expr <- tfcor_expr[!is.na(tfcor_expr$adj.P.Val) & tfcor_expr$adj.P.Val < 0.05 &
                         abs(tfcor_expr$rho) >= 0.4, ]
tfcor_expr <- tfcor_expr[order(-abs(tfcor_expr$rho)), ]
fwrite(tfcor_expr, "results/12_upstream_TF_candidates.csv")
cat(sprintf("[候选] 与 PIEZO1 共表达且 OA 活性差异的 TF: %d 个\n", nrow(tfcor_expr)))

## ============================================================
##  图 fig49 —— TF 活性 OA vs Normal (DoRothEA 差异活性 Top TF)
## ============================================================
topA <- head(ttA[order(ttA$adj.P.Val, -abs(ttA$logFC)), ], 15)
topA <- topA[order(topA$logFC), ]
p49 <- ggplot(topA, aes(x = reorder(tf, logFC), y = logFC, fill = logFC > 0)) +
  geom_bar(stat = "identity", width = 0.7, color = "white") +
  scale_fill_manual(values = c("#2166AC", "#B2182B")) +
  coord_flip() +
  labs(title = "Transcription factor activity: OA vs Normal subchondral bone",
       subtitle = paste0("DoRothEA inference (", ifelse(use_dorothea, "weighted-mean of target expression", "TF expression proxy"), ")"),
       x = NULL, y = "log2(activity OA/Normal)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "none")
ggsave("figs/fig49_TFactivity_OAvsNormal.png", p49, width = 8, height = 6, dpi = 200)

## ============================================================
##  图 fig50 —— 与 PIEZO1 共表达的候选上游 TF (按 |Spearman rho|)
## ============================================================
topTF <- head(cor_all[cor_all$is_TF, ], 20)
topTF <- topTF[order(topTF$rho), ]
p50 <- ggplot(topTF, aes(x = reorder(gene, rho), y = rho, fill = rho > 0)) +
  geom_bar(stat = "identity", width = 0.7, color = "white") +
  scale_fill_manual(values = c("#2166AC", "#B2182B")) +
  coord_flip() +
  labs(title = "Transcription factors co-expressed with PIEZO1",
       subtitle = "Spearman correlation across OA and Normal samples (n=50)",
       x = NULL, y = "Spearman rho (vs PIEZO1 expression)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "none")
ggsave("figs/fig50_PIEZO1_TF_coexpression.png", p50, width = 8, height = 7, dpi = 200)

## ============================================================
##  图 fig51 —— PIEZO1 上游 TF 调控网络
## ============================================================
netTF <- head(cor_all[cor_all$is_TF & abs(cor_all$rho) >= 0.35, ], 15)
if (nrow(netTF) < 8) netTF <- head(cor_all[cor_all$is_TF, ], 12)
edges <- data.frame(from = "PIEZO1", to = netTF$gene, w = abs(netTF$rho))
gnet  <- graph_from_data_frame(edges, directed = FALSE, vertices = c("PIEZO1", netTF$gene))
V(gnet)$type <- ifelse(V(gnet)$name == "PIEZO1", "hub", "tf")
V(gnet)$col  <- ifelse(V(gnet)$name == "PIEZO1", "#000000",
                ifelse(netTF$rho[match(V(gnet)$name, netTF$gene)] > 0, "#B2182B", "#2166AC"))
E(gnet)$w <- edges$w
png("figs/fig51_PIEZO1_TF_network.png", width = 2200, height = 1800, res = 200)
set.seed(42)
lay <- layout_with_fr(gnet, weights = E(gnet)$w)
plot(gnet, layout = lay,
     vertex.size   = c(38, rep(22, vcount(gnet)-1))[order(c(1, match(netTF$gene, V(gnet)$name)))],
     vertex.color  = V(gnet)$col,
     vertex.label  = V(gnet)$name, vertex.label.cex = 1.05,
     vertex.label.color = "white", vertex.label.font = 2,
     vertex.frame.color = "white",
     edge.width    = E(gnet)$w * 10,
     edge.color    = rgb(0.4,0.4,0.4,0.5),
     main = "PIEZO1 upstream transcription-factor network")
legend("topright", legend = c("PIEZO1 (hub)", "co-upregulated TF", "co-downregulated TF"),
       pch = 21, pt.bg = c("#000000","#B2182B","#2166AC"), bty = "n", cex = 1.0)
dev.off()
cat(sprintf("[网络] %d 个 TF 节点接入 PIEZO1\n", nrow(netTF)))

## ============================================================
##  图 fig52 —— PIEZO1 共表达模块 GO 功能注释
## ============================================================
mod_genes <- cor_all$gene[abs(cor_all$rho) >= 0.4]
if (length(mod_genes) < 50) mod_genes <- head(cor_all$gene, 300)
mod_genes <- setdiff(mod_genes, "PIEZO1")
ego <- tryCatch(enrichGO(gene = mod_genes, keyType = "SYMBOL", OrgDb = org.Hs.eg.db,
                         ont = "BP", pAdjustMethod = "BH",
                         pvalueCutoff = 0.05, qvalueCutoff = 0.1),
                error = function(e) NULL)
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  fwrite(as.data.frame(ego), "results/12_PIEZO1_module_GO.csv")
  p52 <- dotplot(ego, showCategory = 15) +
    ggtitle("GO biological processes of the PIEZO1 co-expression module") +
    theme(axis.text.y = element_text(size = 8),
          plot.title = element_text(face = "bold", size = 13))
  ggsave("figs/fig52_PIEZO1_module_GO.png", p52, width = 9, height = 7.5, dpi = 200)
  cat(sprintf("[GO] 共表达模块富集到 %d 条 BP 通路\n", nrow(as.data.frame(ego))))
} else cat("[GO] 共表达模块无显著 BP 通路\n")

## ---------- 汇总 ----------
cat("\n================= 上游 TF 分析汇总 =================\n")
cat(sprintf("TF 活性推断方法 : %s\n", ifelse(use_dorothea, "DoRothEA (TF-target regulons)", "TF 差异表达 (proxy)")))
cat(sprintf("OA vs Normal 差异活性 TF : %d 个 (FDR<0.05)\n", sum(ttA$adj.P.Val < 0.05, na.rm=TRUE)))
cat(sprintf("与 PIEZO1 共表达 TF 候选 : %d 个\n", sum(cor_all$is_TF)))
cat(sprintf("整合上游候选 TF          : %d 个\n", nrow(tfcor_expr)))
cat("输出: figs/fig49~fig52 + results/12_*.csv\n")
cat("===================================================\n")
