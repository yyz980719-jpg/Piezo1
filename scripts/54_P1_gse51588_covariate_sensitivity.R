# =====================================================================
# 54_P1_gse51588_covariate_sensitivity.R
# 模拟审稿后优化 · T6.1 Covariate sensitivity for the GSE51588 models
#
# 审计发现 (重要):
#   1) GSE51588 的 series matrix 实际提供了 age 与 gender —— 原稿与 Supplement
#      均称"协变量不可获得", 必须更正。BMI 与 KL grade 确实未 deposit。
#   2) OA 组平均年龄 69.7 岁, non-OA 组 38.4 岁, 相差约 31 岁;
#      性别分布相近 (55% vs 60% 女性)。年龄是该队列疾病对比的主要潜在混杂。
#   3) donor 编号在 OA 与 non-OA 之间重号 (如 Normal-1 与 OA-1),
#      必须使用 paste(disease, donor) 作为唯一 donor 键 (原脚本 07 已正确处理)。
#
# 本脚本对 PIEZO1 的三个核心对比, 在加 / 不加 age + gender 两种情形下重估。
#   M1 配对 compartment 对比 (OA 内) —— 配对设计自动吸收年龄性别
#   M2 disease x compartment interaction (主要交互检验)
#   M3 disease main effect (region-adjusted)
#
# 输出: results/54_covariate_sensitivity.csv
# =====================================================================

suppressPackageStartupMessages({ library(limma); library(data.table) })
source("scripts/_bootstrap.R")
RES <- PIEZO1_RESULTS
sink(file.path(RES, "54_P1_covariate_sensitivity_log.txt"), split = TRUE)
cat("=== 54_P1 GSE51588 covariate sensitivity ===\n\n")

## ---------------------------------------------------------------------
## 1. 载入表达矩阵与元数据 (含 age / gender)
## ---------------------------------------------------------------------
con <- gzfile("data/GSE51588_series_matrix.txt.gz", "rt")
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

getchar <- function(pattern) {
  ln <- grep("^!Sample_characteristics_ch1", raw, value = TRUE)
  hit <- ln[grepl(pattern, ln, ignore.case = TRUE)]
  if (!length(hit)) return(NULL)
  v <- strsplit(gsub("\"", "", gsub("^!Sample_characteristics_ch1\t", "", hit[1])), "\t")[[1]]
  sub("^[^:]*:\\s*", "", v)
}
age  <- suppressWarnings(as.numeric(getchar("^!Sample_characteristics_ch1\t\"age")))
gend <- getchar("^!Sample_characteristics_ch1\t\"gender")
cat("age available:", sum(!is.na(age)), "/", length(age),
    "| range:", paste(range(age, na.rm = TRUE), collapse = "-"), "\n")
cat("gender available:", sum(!is.na(gend)), "/", length(gend),
    "|", paste(names(table(gend)), table(gend), sep = "=", collapse = " "), "\n")

## 探针注释与基因级合并: 严格复现脚本 07 的做法, 以便与已发表数值对齐
ann <- fread(data_path("GPL13497_probe2symbol.tsv"), header = TRUE, sep = "\t",
             colClasses = "character")
ann <- ann[!is.na(ann$symbol) & nzchar(ann$symbol) &
             grepl("^[A-Za-z][A-Za-z0-9@\\.\\-]*$", ann$symbol)]
common <- intersect(probe, ann$probe)
expr <- expr[probe %in% common, , drop = FALSE]
sym <- ann$symbol[match(probe[probe %in% common], ann$probe)]
## 每个基因保留平均表达量最高的探针 (与脚本 07 一致)
o <- order(-rowMeans(expr, na.rm = TRUE))
expr <- expr[o, , drop = FALSE]; sym <- sym[o]
expr <- expr[o, , drop = FALSE]; sym <- sym[o]
kd <- !duplicated(sym)
expr_l <- expr[kd, , drop = FALSE]; rownames(expr_l) <- sym[kd]
if (max(expr_l, na.rm = TRUE) > 100) expr_l <- log2(expr_l + 1)   # 与脚本 07 相同的尺度规则
cat("genes after collapsing (max-mean probe):", nrow(expr_l), "\n")

parts <- do.call(rbind, strsplit(titles, "-"))
disease <- factor(parts[, 1], levels = c("Normal", "OA"))
comp    <- factor(parts[, 2], levels = c("LT", "MT"))
donor   <- parts[, 3]
uid     <- paste(disease, donor, sep = "_")     # 唯一 donor 键 (OA-1 与 Normal-1 不同人)
meta <- data.table(sample = titles, gsm = colnames(expr), disease = disease,
                   compartment = comp, donor = uid, age = age, gender = factor(gend))
cat("\ndonor keys:", uniqueN(meta$donor), "(若仅用编号则为", uniqueN(donor), ", 会错误合并 OA 与 non-OA 供体)\n")
cat("\nmetadata by disease group:\n")
print(meta[, .(n = .N, n_donors = uniqueN(donor),
               age_mean = round(mean(age, na.rm = TRUE), 1),
               age_sd = round(sd(age, na.rm = TRUE), 1),
               pct_female = round(100 * mean(gender == "Female", na.rm = TRUE), 1)), by = disease])
tt_age <- t.test(age ~ disease, data = meta)
cat(sprintf("\n年龄差异 OA - Normal = %+.1f 岁 (Welch t=%.2f, P=%.3g)\n",
            diff(tapply(meta$age, meta$disease, mean))[["OA"]], tt_age$statistic, tt_age$p.value))
cat(">> OA 与 non-OA 组年龄严重不匹配; 这是 disease 相关推断的主要混杂来源。\n\n")

y <- expr_l["PIEZO1", ]
out <- list()
add <- function(model, contrast, adjustment, est, stat, p, note = "") {
  out[[length(out) + 1]] <<- data.table(model = model, contrast = contrast,
    adjustment = adjustment, estimate = est, statistic = stat, P = p, note = note)
}

## ---------------------------------------------------------------------
## M1 配对 compartment 对比 (OA 内) —— 配对设计自动吸收 age/gender
## ---------------------------------------------------------------------
oa <- meta[disease == "OA"]
ids <- intersect(oa$donor[oa$compartment == "LT"], oa$donor[oa$compartment == "MT"])
## 必须按 donor 排序后再配对, 否则 t.test 会按输入顺序错配
da <- oa[oa$compartment == "LT" & oa$donor %in% ids][order(donor)]
db <- oa[oa$compartment == "MT" & oa$donor %in% ids][order(donor)]
stopifnot(identical(da$donor, db$donor))
a <- y[match(da$gsm, colnames(expr))]
b <- y[match(db$gsm, colnames(expr))]
tt <- t.test(a, b, paired = TRUE)
add("M1 paired compartment (within OA)", "lateral - medial", "none (paired within donor)",
    tt$estimate, tt$statistic, tt$p.value,
    "donor-invariant covariates (age, gender) are absorbed by the pairing")

## ---------------------------------------------------------------------
## M2 disease x compartment interaction (donor 层面的 medial-lateral delta)
## ---------------------------------------------------------------------
uids <- unique(meta$donor)
dl <- lapply(uids, function(u) {
  m <- meta$compartment[meta$donor == u] == "MT"
  l <- meta$compartment[meta$donor == u] == "LT"
  if (sum(m) == 1 && sum(l) == 1)
    y[match(meta$gsm[meta$donor == u][m], colnames(expr))] -
    y[match(meta$gsm[meta$donor == u][l], colnames(expr))]
  else NULL
})
names(dl) <- uids
dl <- dl[!sapply(dl, is.null)]
D  <- do.call(cbind, dl)
rownames(D) <- "PIEZO1"     # cbind 会沿用 GSM 名, 显式改回基因名
dinfo <- data.table(donor = names(dl),
                    disease = factor(sub("_.*", "", names(dl)), levels = c("Normal", "OA")),
                    age = meta$age[match(names(dl), meta$donor)],
                    gender = meta$gender[match(names(dl), meta$donor)])
cat("donor-level deltas:", ncol(D), "| OA:", sum(dinfo$disease == "OA"),
    "| Normal:", sum(dinfo$disease == "Normal"), "\n")

## (a) 原稿主检验: moderated limma, 无协变量 (P = 0.147)
d0 <- model.matrix(~ disease, data = dinfo); colnames(d0) <- c("(Intercept)", "disOA")
f0 <- eBayes(lmFit(D, d0), trend = TRUE)
add("M2 disease x compartment interaction", "delta(MT-LT): OA - Normal",
    "limma moderated, no covariate", f0$coefficients[1, "disOA"],
    f0$t[1, "disOA"], f0$p.value[1, "disOA"],
    "primary interaction test reported in the manuscript (P=0.147)")

## (b) 协变量敏感性: moderated limma + age + gender
d1 <- model.matrix(~ disease + age + gender, data = dinfo)
colnames(d1)[colnames(d1) == "diseaseOA"] <- "disOA"
f1 <- eBayes(lmFit(D, d1), trend = TRUE)
add("M2 disease x compartment interaction", "delta(MT-LT): OA - Normal",
    "limma moderated + age + gender", f1$coefficients[1, "disOA"],
    f1$t[1, "disOA"], f1$p.value[1, "disOA"],
    "covariate-adjusted sensitivity (age and gender ARE deposited)")

## (c) 非配对 t 检验 (原 P = 0.064) 与 + age + gender
t0 <- t.test(as.numeric(D[1, ]) ~ dinfo$disease)
add("M2 unmoderated sensitivity", "delta(MT-LT): OA - Normal", "none (Welch/Student t)",
    t0$estimate[2] - t0$estimate[1], t0$statistic, t0$p.value,
    "unmoderated sensitivity (P=0.064 in the manuscript)")
g1 <- lm(as.numeric(D[1, ]) ~ disease + age + gender, data = dinfo)
s1 <- summary(g1)$coefficients["diseaseOA", ]
add("M2 unmoderated sensitivity", "delta(MT-LT): OA - Normal", "+ age + gender",
    s1[1], s1[3], s1[4], "covariate-adjusted sensitivity")

## ---------------------------------------------------------------------
## M3 disease main effect (region-adjusted)
## ---------------------------------------------------------------------
l0 <- lm(y ~ disease + compartment, data = meta)
s2 <- summary(l0)$coefficients["diseaseOA", ]
add("M3 disease main effect", "OA - Normal (region adjusted)", "compartment",
    s2[1], s2[3], s2[4], "null result")
l1 <- lm(y ~ disease + compartment + age + gender, data = meta)
s3 <- summary(l1)$coefficients["diseaseOA", ]
add("M3 disease main effect", "OA - Normal (region adjusted)", "compartment + age + gender",
    s3[1], s3[3], s3[4], "covariate-adjusted sensitivity")

res <- rbindlist(out)
cat("\n=== 结果 ===\n")
print(res[, .(model, contrast, adjustment, estimate = round(estimate, 4),
              stat = round(statistic, 3), P = signif(P, 4))], nrows = 30)
cat("\n说明: M1 为配对设计, age/gender 在 donor 内不变, 配对本身即完成控制。\n")
cat("      M2/M3 为 donor 间比较, 已用 GSE51588 提供的 age 与 gender 做校正敏感性分析。\n")
cat("      BMI 与 KL grade 未被 deposit, 无法校正, 列为 limitation。\n")
fwrite(res, file.path(RES, "54_covariate_sensitivity.csv"))
cat("\n输出: results/54_covariate_sensitivity.csv\n")
sink()
