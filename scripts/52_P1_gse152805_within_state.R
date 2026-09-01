# =====================================================================
# 52_P1_gse152805_within_state.R
# 模拟审稿后优化 · T3.1/T3.2 软骨内 PIEZO1 方向不一致的解析
#
# 问题: GSE104782 (单细胞) 显示 PIEZO1 在 homeostatic 端更高;
#       但 GSE57218 (preserved-affected) 与 GSE152805 (lateral-medial)
#       都把更高的 PIEZO1 放在 affected/medial 一侧。
#       审稿人要求区分: 总体差异来自 (a) 细胞状态构成改变, 还是
#       (b) 同一状态内部的 PIEZO1 表达改变。
#
# 设计: 在 GSE152805 的 3 个配对 donor 内, 按 chondrocyte state 分层,
#       计算每个 state 内 lateral - medial 的 PIEZO1 差异 (donor 层面报告方向,
#       不用 cell-level P 值做主证据); 并做 Kitagawa 分解, 把总体差异拆成
#       within-state 成分与 composition 成分。
#
# 输出: results/52_*.csv
# =====================================================================

suppressPackageStartupMessages({ library(Matrix); library(data.table); library(Seurat) })
source("scripts/_bootstrap.R")
RES <- PIEZO1_RESULTS
set.seed(1)
sink(file.path(RES, "52_P1_within_state_log.txt"), split = TRUE)
cat("=== 52_P1 GSE152805 within-state PIEZO1 ===\n\n")

## ---------------------------------------------------------------------
## 1. 载入 GSE152805 并 QC (与脚本 39 完全一致)
## ---------------------------------------------------------------------
dir <- data_path("GSE152805_raw")
files <- list.files(dir)
keep <- c("OA_oLT_113","OA_oLT_116","OA_oLT_118","OA_MT_113","OA_MT_116","OA_MT_118")
counts <- NULL; meta <- data.frame()
for (s in keep) {
  f_m <- file.path(dir, list.files(dir, pattern = paste0("^GSM[0-9]+_", s, "\\.matrix")))
  f_b <- file.path(dir, list.files(dir, pattern = paste0("^GSM[0-9]+_", s, "\\.barcodes")))
  f_g <- file.path(dir, list.files(dir, pattern = paste0("^GSM[0-9]+_", s, "\\.genes")))
  M <- readMM(gzfile(f_m))
  bc <- read.delim(gzfile(f_b), header = FALSE, stringsAsFactors = FALSE)[, 1]
  gn <- read.delim(gzfile(f_g), header = FALSE, stringsAsFactors = FALSE)
  dimnames(M) <- list(gn[, 2], bc)
  cat(sprintf("  %s: %d genes x %d cells\n", s, nrow(M), ncol(M)))
  counts <- if (is.null(counts)) M else cbind(counts, M)
  comp <- if (grepl("oLT", s)) "lateral_nonload" else "medial_load"
  don  <- sub(".*_", "", s)
  meta <- rbind(meta, data.frame(cell = paste(s, bc, sep = "_"),
                                 compartment = comp, donor = don, stringsAsFactors = FALSE))
}
rownames(meta) <- meta$cell
meta$compartment <- factor(meta$compartment, levels = c('medial_load','lateral_nonload'))
cat("combined:", dim(counts), "\n")

nc <- Matrix::colSums(counts); nf <- Matrix::colSums(counts > 0)
mt <- grep("^MT-", rownames(counts), value = TRUE)
mp <- if (length(mt)) Matrix::colSums(counts[mt, , drop = FALSE]) / nc else rep(0, ncol(counts))
ok <- nf >= 200 & nc >= 500 & mp < 0.25
cat("cells passing QC:", sum(ok), "of", ncol(counts), "\n")
counts <- counts[, ok, drop = FALSE]; meta <- meta[ok, ]
if (any(duplicated(rownames(counts)))) counts <- counts[!duplicated(rownames(counts)), ]
lib  <- Matrix::colSums(counts)
expr <- Matrix::t(log1p(Matrix::t(counts) / lib * 1e4))
expr <- as.matrix(expr)
cat("expression matrix:", dim(expr), "\n\n")

## ---------------------------------------------------------------------
## 2. Programme signature 与细胞状态分配 (最大投影得分 = 最相似状态)
## ---------------------------------------------------------------------
seu <- readRDS("data/scrna/GSE104782_seurat.rds")
sub <- as.character(seu@meta.data$subtype)
d0  <- tryCatch(GetAssayData(seu, layer = "data"), error = function(e) GetAssayData(seu, slot = "data"))
lin <- as.matrix(expm1(d0)); subtypes <- sort(unique(sub))
mean_mat <- sapply(subtypes, function(s) rowMeans(lin[, sub == s, drop = FALSE]))
sig_genes <- unique(unlist(lapply(subtypes, function(s) {
  o <- rowMeans(mean_mat[, setdiff(subtypes, s), drop = FALSE])
  fc <- (mean_mat[, s] + 1e-6) / (o + 1e-6); fc[mean_mat[, s] < 0.05] <- 0
  names(sort(fc, decreasing = TRUE))[1:100]
})))
sig <- mean_mat[sig_genes, , drop = FALSE]
sig_z <- t(scale(t(sig))); sig_z[is.na(sig_z)] <- 0
common <- intersect(sig_genes, rownames(expr))
zz <- t(scale(t(expr[common, , drop = FALSE]))); zz[is.na(zz)] <- 0
W <- t(zz) %*% sig_z[common, , drop = FALSE]
colnames(W) <- subtypes
cat("cells scored:", nrow(W), "\n")

state <- colnames(W)[max.col(W, ties.method = "first")]
state[max.col(W, ties.method = "first") == 0] <- NA
## 低置信分配: 最高分 <= 0 视为不可分配
state[apply(W, 1, max) <= 0] <- "unassigned"
meta$state <- state
pz <- if ("PIEZO1" %in% rownames(expr)) as.numeric(expr["PIEZO1", ]) else rep(NA_real_, ncol(expr))
meta$PIEZO1 <- pz
cat("\nstate assignment counts:\n"); print(table(meta$state, meta$compartment))
cat("\nPIEZO1 detection overall:", round(mean(pz > 0) * 100, 1), "%\n\n")

## ---------------------------------------------------------------------
## 3. 每个 donor x state 的 lateral - medial 差异
## ---------------------------------------------------------------------
MIN_CELLS <- 50     # 每个 donor x compartment x state 的最小细胞数
dons <- sort(unique(meta$donor))
rows <- list()
for (dd in dons) for (s in subtypes) {
  a <- meta$PIEZO1[meta$donor == dd & meta$state == s & meta$compartment == "lateral_nonload"]
  b <- meta$PIEZO1[meta$donor == dd & meta$state == s & meta$compartment == "medial_load"]
  na <- length(a); nb <- length(b)
  if (na < MIN_CELLS || nb < MIN_CELLS) {
    rows[[length(rows) + 1]] <- data.table(donor = dd, state = s, n_lateral = na, n_medial = nb,
      mean_lateral = NA_real_, mean_medial = NA_real_, delta_lateral_minus_medial = NA_real_,
      det_lateral = NA_real_, det_medial = NA_real_, available = FALSE)
  } else {
    rows[[length(rows) + 1]] <- data.table(donor = dd, state = s, n_lateral = na, n_medial = nb,
      mean_lateral = mean(a), mean_medial = mean(b),
      delta_lateral_minus_medial = mean(a) - mean(b),
      det_lateral = mean(a > 0), det_medial = mean(b > 0), available = TRUE)
  }
}
within <- rbindlist(rows)
cat("=== within-state PIEZO1: lateral - medial, per donor ===\n")
print(within[, .(donor, state, n_lateral, n_medial,
                 mean_lateral = round(mean_lateral, 4), mean_medial = round(mean_medial, 4),
                 delta = round(delta_lateral_minus_medial, 4),
                 available)], nrows = 30)

av <- within[available == TRUE]
cat("\n=== 方向一致性 (仅 available 的 donor x state 组合) ===\n")
summ <- av[, .(n_donors = .N, n_negative = sum(delta_lateral_minus_medial < 0),
               mean_delta = mean(delta_lateral_minus_medial)), by = state]
summ[, direction := ifelse(n_negative == n_donors, "medial higher in all donors",
                    ifelse(n_negative == 0, "lateral higher in all donors", "inconsistent"))]
print(summ[, .(state, n_donors, n_negative, mean_delta = round(mean_delta, 4), direction)])

## ---------------------------------------------------------------------
## 4. Kitagawa 分解: 总体 PIEZO1 差异 = within-state 成分 + composition 成分
##    delta_total = sum_s( w_bar_s * delta_within_s ) + sum_s( delta_w_s * mu_bar_s )
##    w_bar_s = 两侧状态占比均值; delta_w_s = w_s(lateral) - w_s(medial)
## ---------------------------------------------------------------------
dec <- list()
for (dd in dons) {
  m <- meta[meta$donor == dd, ]
  la <- m[m$compartment == "lateral_nonload", ]; me <- m[m$compartment == "medial_load", ]
  wl <- table(factor(la$state, levels = subtypes)) / nrow(la)
  wm <- table(factor(me$state, levels = subtypes)) / nrow(me)
  mu_l <- sapply(subtypes, function(s) mean(la$PIEZO1[la$state == s]))
  mu_m <- sapply(subtypes, function(s) mean(me$PIEZO1[me$state == s]))
  ok <- !is.na(mu_l) & !is.na(mu_m)
  wbar <- (as.numeric(wl) + as.numeric(wm)) / 2
  dw   <- as.numeric(wl) - as.numeric(wm)
  dwithin <- mu_l - mu_m
  mubar <- (mu_l + mu_m) / 2
  within_part <- sum(wbar[ok] * dwithin[ok], na.rm = TRUE)
  comp_part   <- sum(dw[ok] * mubar[ok], na.rm = TRUE)
  total <- mean(la$PIEZO1) - mean(me$PIEZO1)
  dec[[length(dec) + 1]] <- data.table(donor = dd,
    total_delta_lateral_minus_medial = total,
    within_state_component = within_part,
    composition_component = comp_part,
    sum_check = within_part + comp_part)
}
dec <- rbindlist(dec)
dec[, pct_within := 100 * within_state_component / (within_state_component + composition_component)]
cat("\n=== Kitagawa 分解 (总体差异 = within-state + composition) ===\n")
print(dec[, .(donor, total = round(total_delta_lateral_minus_medial, 4),
              within = round(within_state_component, 4), comp = round(composition_component, 4),
              check = round(sum_check, 4), pct_within = round(pct_within, 1))])
cat("\nmean within-state contribution: ",
    round(mean(dec$within_state_component), 4),
    " | mean composition contribution: ", round(mean(dec$composition_component), 4), "\n")

## ---------------------------------------------------------------------
## 5. 总体 (未分层) 的 donor-level 方向 —— 与 within-state 结果对照
## ---------------------------------------------------------------------
ov <- rbindlist(lapply(dons, function(dd) {
  a <- meta$PIEZO1[meta$donor == dd & meta$compartment == "lateral_nonload"]
  b <- meta$PIEZO1[meta$donor == dd & meta$compartment == "medial_load"]
  data.table(donor = dd, n_lateral = length(a), n_medial = length(b),
             overall_delta_lateral_minus_medial = mean(a) - mean(b),
             det_lateral = mean(a > 0), det_medial = mean(b > 0))
}))
cat("\n=== 总体 (未分层) 方向 ===\n"); print(ov, digits = 4)
cat("donors with medial PIEZO1 higher overall:", sum(ov$overall_delta_lateral_minus_medial < 0),
    "/", nrow(ov), "\n")

## 状态构成 (lateral - medial, 百分点)
comp <- rbindlist(lapply(dons, function(dd) {
  m <- meta[meta$donor == dd, ]
  wl <- table(factor(m$state[m$compartment == "lateral_nonload"], levels = subtypes))
  wm <- table(factor(m$state[m$compartment == "medial_load"], levels = subtypes))
  wl <- wl / sum(wl); wm <- wm / sum(wm)
  data.table(donor = dd, state = subtypes,
             prop_lateral = as.numeric(wl), prop_medial = as.numeric(wm),
             delta_prop = as.numeric(wl) - as.numeric(wm))
}))
cat("\n=== 状态构成: lateral - medial (比例差) ===\n")
print(dcast(comp, state ~ donor, value.var = "delta_prop"), digits = 3)

fwrite(within, file.path(RES, "52_within_state_per_donor.csv"))
fwrite(dec,    file.path(RES, "52_kitagawa_decomposition.csv"))
fwrite(ov,     file.path(RES, "52_overall_direction.csv"))
fwrite(comp,   file.path(RES, "52_state_composition.csv"))
cat("\n输出: results/52_*.csv\n")
sink()
