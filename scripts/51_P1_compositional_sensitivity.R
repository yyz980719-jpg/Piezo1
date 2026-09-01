# =====================================================================
# 51_P1_compositional_sensitivity.R
# 模拟审稿后优化 · T2.1 Programme weights 的组成型 (compositional) 敏感性
#
# 背景: NNLS 拟合得到的 7 个 programme weights 被重新缩放为 sum = 1。
#       闭合数据 (closed/compositional data) 存在"一个分量上升必然压低其它分量"
#       的数学耦合。审稿人要求证明 HomC↑ / preHTC↓ 不只是闭合约束的产物。
#
# 策略 (5 个层次, 层层剥离闭合约束):
#   L0 原始 sum-to-one 权重                     —— 现稿做法 (闭合)
#   L1 CLR (中心化 log 比, 多种 pseudocount)     —— 闭合一致变换
#   L2 ALR (以 7 个分量逐一为参照)                —— 闭合一致变换
#   L3 log(HomC / preHTC) 单一统计量             —— 对闭合与缩放完全不变 (主统计量)
#   L4 unclosed gene-set score (各 programme 单独打分, 不跨 programme 归一)
#                                               —— 完全无闭合约束 (第二独立证据)
#
# 适用数据集: GSE51588 (软骨下骨, 配对 medial/lateral)
#             GSE57218 (关节软骨, 配对 preserved/OA-affected)
# 注: GSE152805 使用的是无约束投影分数 (W = z(expr) %*% z(sig)), 本身不受闭合约束,
#     因此不在此脚本内重复分析, 仅在稿件中明确其量纲与 bulk 权重不同。
#
# 输出: results/51_*.csv
# =====================================================================

suppressPackageStartupMessages({ library(Matrix); library(data.table); library(Seurat) })
source("scripts/_bootstrap.R")
RES <- PIEZO1_RESULTS
set.seed(1)
sink(result_path("51_P1_compositional_sensitivity_log.txt"), split = TRUE)
cat("=== 51_P1 programme compositional sensitivity ===\n\n")

## ---------------------------------------------------------------------
## 0. 从 GSE104782 重建 programme signature (与 33/35 完全一致)
## ---------------------------------------------------------------------
seu <- readRDS(data_path("scrna", "GSE104782_seurat.rds"))
sub <- as.character(seu@meta.data$subtype)
d0  <- tryCatch(Seurat::GetAssayData(seu, layer = "data"),
                error = function(e) Seurat::GetAssayData(seu, slot = "data"))
lin <- as.matrix(expm1(d0)); subtypes <- sort(unique(sub))
cat("GSE104782 cells:", ncol(seu), "| subtypes:", paste(subtypes, collapse = "/"), "\n")

mean_mat <- sapply(subtypes, function(s) rowMeans(lin[, sub == s, drop = FALSE]))
sig_list <- lapply(subtypes, function(s) {
  others <- rowMeans(mean_mat[, setdiff(subtypes, s), drop = FALSE])
  fc <- (mean_mat[, s] + 1e-6) / (others + 1e-6)
  fc[mean_mat[, s] < 0.05] <- 0
  names(sort(fc, decreasing = TRUE))[1:100]
})
names(sig_list) <- subtypes
sig_genes <- unique(unlist(sig_list))
sig <- mean_mat[sig_genes, , drop = FALSE]
sig_z <- t(scale(t(sig))); sig_z[is.na(sig_z)] <- 0
cat("signature genes (union):", length(sig_genes),
    "| per programme:", length(sig_list[[1]]), "\n\n")

## ---------------------------------------------------------------------
## 1. 组成型变换工具
## ---------------------------------------------------------------------
## multiplicative simple replacement: 零值用 0.5*最小正值替换, 并从非零分量等比扣除
msr <- function(x, delta = NULL) {
  x <- as.numeric(x); z <- x <= 0
  if (all(z)) return(rep(1/length(x), length(x)))
  if (is.null(delta)) delta <- 0.5 * min(x[!z])
  if (!any(z)) return(x)
  xr <- x; xr[z] <- delta
  s <- sum(x[z])                       # 需要扣除的总量
  xr[!z] <- xr[!z] * (1 - s / sum(xr[!z]))
  xr
}
clr <- function(x) { lx <- log(x); lx - mean(lx) }
alr <- function(x, ref) { lx <- log(x); lx - lx[ref] }

## NNLS (与 33/35 相同的实现)
nnls_solve <- function(A, b) {
  r <- optim(rep(1/ncol(A), ncol(A)),
             function(x) sum((as.vector(A %*% x) - b)^2),
             function(x) as.vector(2 * t(A) %*% (as.vector(A %*% x) - b)),
             method = "L-BFGS-B", lower = rep(0, ncol(A)),
             control = list(maxit = 3000, factr = 1e7))
  x <- r$par; if (sum(x) > 0) x <- x / sum(x); x
}

## ---------------------------------------------------------------------
## 2. 通用: 给定 per-sample 权重矩阵与配对信息, 做 5 层敏感性分析
## ---------------------------------------------------------------------
layered_test <- function(W, grp_a, grp_b, donor, label, score_mat = NULL) {
  ## grp_a / grp_b: 每个样本所属组别; donor: 配对 id; 统计量 = mean(a) - mean(b)
  ids <- intersect(donor[grp_a], donor[grp_b])
  ia <- which(grp_a & donor %in% ids); ib <- which(grp_b & donor %in% ids)
  ia <- ia[order(donor[ia])]; ib <- ib[order(donor[ib])]
  stopifnot(length(ia) > 1L, length(ia) == length(ib), all(donor[ia] == donor[ib]))

  ptest <- function(v) {
    tt <- t.test(v[ia], v[ib], paired = TRUE)
    data.frame(delta = mean(v[ia] - v[ib]), t = tt$statistic, p = tt$p.value)
  }
  out <- list()

  ## L0 原始 sum-to-one 权重
  l0 <- do.call(rbind, lapply(subtypes, function(s) {
    r <- ptest(W[, s]); data.frame(dataset = label, layer = "L0_raw_sum1",
      programme = s, delta = r$delta, t = r$t, p = r$p) }))
  out[[length(out) + 1]] <- l0

  ## L1 CLR (三种 zero 处理)
  for (zd in c("msr_0.5min", "add_1e-3", "add_1e-2")) {
    Wr <- t(apply(W, 1, function(row) switch(zd,
      "msr_0.5min" = msr(row),
      "add_1e-3"   = row + 1e-3,
      "add_1e-2"   = row + 1e-2)))
    colnames(Wr) <- subtypes
    Wr <- Wr / rowSums(Wr)
    Wc <- t(apply(Wr, 1, clr)); colnames(Wc) <- subtypes
    out[[length(out) + 1]] <- do.call(rbind, lapply(subtypes, function(s) {
      r <- ptest(Wc[, s]); data.frame(dataset = label, layer = paste0("L1_CLR_", zd),
        programme = s, delta = r$delta, t = r$t, p = r$p) }))
  }

  ## L2 ALR (7 个参照逐一)
  Wrm <- t(apply(W, 1, msr)); colnames(Wrm) <- subtypes
  Wrm <- Wrm / rowSums(Wrm)
  for (ref in subtypes) {
    Wa <- t(apply(Wrm, 1, alr, ref = ref)); colnames(Wa) <- subtypes
    out[[length(out) + 1]] <- do.call(rbind, lapply(setdiff(subtypes, ref), function(s) {
      r <- ptest(Wa[, s]); data.frame(dataset = label, layer = paste0("L2_ALR_ref_", ref),
        programme = s, delta = r$delta, t = r$t, p = r$p) }))
  }

  ## L3 关键闭合不变量: log(HomC / preHTC)
  lr <- log(Wrm[, "HomC"]) - log(Wrm[, "preHTC"])
  r <- ptest(lr)
  out[[length(out) + 1]] <- data.frame(dataset = label, layer = "L3_logratio_HomC_over_preHTC",
    programme = "log(HomC/preHTC)", delta = r$delta, t = r$t, p = r$p)

  ## L4 unclosed gene-set score (无闭合约束)
  if (!is.null(score_mat)) {
    l4 <- do.call(rbind, lapply(subtypes, function(s) {
      r <- ptest(score_mat[, s]); data.frame(dataset = label, layer = "L4_unclosed_gs_score",
        programme = s, delta = r$delta, t = r$t, p = r$p) }))
    out[[length(out) + 1]] <- l4
    r2 <- ptest(score_mat[, "HomC"] - score_mat[, "preHTC"])
    out[[length(out) + 1]] <- data.frame(dataset = label,
      layer = "L4_unclosed_gs_score", programme = "HomC_minus_preHTC",
      delta = r2$delta, t = r2$t, p = r2$p)
  }
  rbindlist(out)
}

## ---------------------------------------------------------------------
## 3. GSE51588 (配对 medial(MT) vs lateral(LT))
## ---------------------------------------------------------------------
cat("--- GSE51588 ---\n")
con <- gzfile(data_path("GSE51588_series_matrix.txt.gz"), "rt")
raw <- readLines(con); close(con)
tstart <- grep("^!series_matrix_table_begin", raw); tend <- grep("^!series_matrix_table_end", raw)
titles <- strsplit(gsub("\"", "", gsub("^!Sample_title\t", "", raw[grep("^!Sample_title", raw)])), "\t")[[1]]
body <- raw[(tstart + 2):(tend - 1)]; body <- body[!grepl("^!", body)]
hdr <- strsplit(gsub("\"", "", raw[tstart + 1]), "\t")[[1]]
tab <- read.delim(text = paste(body, collapse = "\n"), header = FALSE, quote = "",
                  stringsAsFactors = FALSE, check.names = FALSE)
colnames(tab) <- hdr
idcol <- grep("^ID_REF$|^ID_ref$", hdr, ignore.case = TRUE)[1]
probe <- gsub("\"", "", as.character(tab[[idcol]]))
expr <- as.matrix(tab[, setdiff(seq_along(hdr), idcol), drop = FALSE])
expr <- matrix(suppressWarnings(as.numeric(gsub("\"", "", expr))), nrow = nrow(expr))
colnames(expr) <- hdr[setdiff(seq_along(hdr), idcol)]
p2s <- read.delim(data_path("GPL13497_probe2symbol.tsv"), stringsAsFactors = FALSE)
sym <- p2s$symbol[match(probe, p2s$probe)]
keep <- !is.na(sym) & sym != ""; expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
tb <- table(sym)
expr_g <- rowsum(expr, group = sym, reorder = TRUE) / as.numeric(tb[rownames(rowsum(expr, group = sym, reorder = TRUE))])
if (max(expr_g, na.rm = TRUE) < 60) expr_g <- 2^expr_g
parts <- do.call(rbind, strsplit(titles, "-"))
disease <- parts[, 1]; comp <- parts[, 2]; donor <- parts[, 3]
cat("samples:", ncol(expr_g), "| compartments:", paste(unique(comp), collapse = "/"), "\n")

common <- intersect(rownames(sig_z), rownames(expr_g)); S <- sig_z[common, , drop = FALSE]
Zall <- t(scale(t(expr_g[common, ]))); Zall[is.na(Zall)] <- 0
W515 <- t(sapply(seq_len(ncol(expr_g)), function(j) nnls_solve(S, Zall[, j])))
colnames(W515) <- subtypes
cat("zero cells in weights:", sum(W515 == 0), "of", length(W515), "\n")

## L4: 各 programme 独立 gene-set score (mean z across its own top-100 genes)
gs515 <- sapply(subtypes, function(s) {
  g <- intersect(sig_list[[s]], rownames(Zall))
  colMeans(Zall[g, , drop = FALSE], na.rm = TRUE)
})
res515 <- layered_test(W515, comp == "LT", comp == "MT", donor,
                       "GSE51588(LT-MT)", score_mat = gs515)

## ---------------------------------------------------------------------
## 4. GSE57218 (配对 preserved vs OA-affected)
## ---------------------------------------------------------------------
cat("\n--- GSE57218 ---\n")
con <- gzfile(data_path("GSE57218_series_matrix.txt.gz"), "rt")
raw <- readLines(con); close(con)
tstart <- grep("^!series_matrix_table_begin", raw); tend <- grep("^!series_matrix_table_end", raw)
titles <- strsplit(gsub("\"", "", gsub("^!Sample_title\t", "", raw[grep("^!Sample_title", raw)])), "\t")[[1]]
body <- raw[(tstart + 2):(tend - 1)]; body <- body[!grepl("^!", body)]
hdr <- strsplit(gsub("\"", "", raw[tstart + 1]), "\t")[[1]]
tab <- read.delim(text = paste(body, collapse = "\n"), header = FALSE, quote = "",
                  stringsAsFactors = FALSE, check.names = FALSE)
colnames(tab) <- hdr
idcol <- grep("^ID_REF$|^ID$", hdr, ignore.case = TRUE)[1]
probe <- gsub("\"", "", as.character(tab[[idcol]]))
expr <- as.matrix(tab[, setdiff(seq_along(hdr), idcol), drop = FALSE])
expr <- matrix(suppressWarnings(as.numeric(gsub("\"", "", expr))), nrow = nrow(expr))
colnames(expr) <- hdr[setdiff(seq_along(hdr), idcol)]
zf <- gzfile(data_path("GPL6947.annot.gz"), "rt"); an <- readLines(zf); close(zf)
an <- an[!grepl("^\\^|^!", an)]
ah <- which(grepl("^ID\t", an))[1]
ahdr <- strsplit(an[ah], "\t")[[1]]
i_id <- grep("^ID$|^ID_REF$|^Probe", ahdr)[1]
i_sym <- grep("^Gene symbol$|^Symbol$", ahdr)[1]
sp <- strsplit(an[(ah + 1):length(an)], "\t", fixed = TRUE)
pid <- sapply(sp, function(x) if (length(x) >= max(i_id, i_sym)) trimws(x[i_id]) else NA)
psy <- sapply(sp, function(x) if (length(x) >= max(i_id, i_sym)) trimws(x[i_sym]) else NA)
p2s <- data.frame(probe = pid, symbol = psy, stringsAsFactors = FALSE)
p2s <- p2s[!is.na(p2s$symbol) & p2s$symbol != "" & p2s$symbol != "NA", ]
sym <- p2s$symbol[match(probe, p2s$probe)]
keep <- !is.na(sym) & sym != ""; expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
tb <- table(sym)
expr_g2 <- rowsum(expr, group = sym, reorder = TRUE) / as.numeric(tb[rownames(rowsum(expr, group = sym, reorder = TRUE))])
raak <- sub("cartilage_RAAK_", "", titles); raak <- sub("_.*$", "", raak)
ttype <- sub(".*_", "", titles)
cat("samples:", ncol(expr_g2), "| types:", paste(names(table(ttype)), table(ttype), sep = "=", collapse = " "), "\n")

common2 <- intersect(rownames(sig_z), rownames(expr_g2)); S2 <- sig_z[common2, , drop = FALSE]
Z2 <- t(scale(t(expr_g2[common2, ]))); Z2[is.na(Z2)] <- 0
W572 <- t(sapply(seq_len(ncol(expr_g2)), function(j) nnls_solve(S2, Z2[, j])))
colnames(W572) <- subtypes
cat("zero cells in weights:", sum(W572 == 0), "of", length(W572), "\n")
gs572 <- sapply(subtypes, function(s) {
  g <- intersect(sig_list[[s]], rownames(Z2))
  colMeans(Z2[g, , drop = FALSE], na.rm = TRUE)
})
res572 <- layered_test(W572, ttype == "Preserved", ttype == "OA", raak,
                       "GSE57218(Pres-OA)", score_mat = gs572)

## ---------------------------------------------------------------------
## 5. 汇总
## ---------------------------------------------------------------------
all <- rbind(res515, res572)
all[, p_adj_within_layer := p.adjust(p, "BH"), by = .(dataset, layer)]

cat("\n\n=== 关键结果: HomC 与 preHTC 在 5 层分析下的表现 ===\n")
key <- all[programme %in% c("HomC", "preHTC", "log(HomC/preHTC)", "HomC_minus_preHTC")]
key <- key[layer %in% c("L0_raw_sum1", "L1_CLR_msr_0.5min", "L1_CLR_add_1e-2",
                        "L3_logratio_HomC_over_preHTC", "L4_unclosed_gs_score")]
print(key[order(dataset, programme, layer),
          .(dataset, layer, programme, delta = round(delta, 4), p = signif(p, 3))],
      nrows = 60)

cat("\n=== ALR 稳健性: HomC 在 7 个不同参照下的 delta 与 P ===\n")
alr_h <- all[grepl("^L2_ALR", layer) & programme == "HomC"]
print(alr_h[, .(dataset, layer, delta = round(delta, 4), p = signif(p, 4))], nrows = 30)
cat("\nALR HomC P 值范围: ",
    paste(alr_h[, .(min = signif(min(p), 3), max = signif(max(p), 3)), by = dataset][,
               paste0(dataset, " ", min, "-", max)], collapse = " | "), "\n")

fwrite(all, file.path(RES, "51_compositional_sensitivity_full.csv"))
fwrite(all[programme %in% c("HomC", "preHTC", "log(HomC/preHTC)", "HomC_minus_preHTC")],
       file.path(RES, "51_compositional_sensitivity_key.csv"))
fwrite(as.data.table(W515), file.path(RES, "51_weights_GSE51588.csv"))
fwrite(as.data.table(W572), file.path(RES, "51_weights_GSE57218.csv"))
cat("\n输出: results/51_*.csv\n")
sink()
