# 33: scRNA-informed deconvolution of paired bulk + mechanotransduction scoring
source("scripts/_bootstrap.R")
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
set.seed(1); RES <- PIEZO1_RESULTS

cat("== scRNA reference (GSE104782) ==\n")
seu <- readRDS("data/scrna/GSE104782_seurat.rds")
md  <- seu@meta.data; sub <- as.character(md$subtype)
cat("cells:", ncol(seu), "genes:", nrow(seu), "\n"); print(table(sub))
dat <- tryCatch(GetAssayData(seu, layer = "data"), error = function(e) GetAssayData(seu, slot = "data"))
lin <- as.matrix(expm1(dat))
subtypes <- sort(unique(sub))
mean_mat <- sapply(subtypes, function(s) rowMeans(lin[, sub == s, drop = FALSE]))
sig_genes <- unique(unlist(lapply(subtypes, function(s) {
  others <- rowMeans(mean_mat[, setdiff(subtypes, s), drop = FALSE])
  fc <- (mean_mat[, s] + 1e-6) / (others + 1e-6)
  fc[mean_mat[, s] < 0.05] <- 0
  names(sort(fc, decreasing = TRUE))[1:100]
})))
cat("signature genes:", length(sig_genes), "\n")
sig <- mean_mat[sig_genes, , drop = FALSE]
sig_z <- t(scale(t(sig))); sig_z[is.na(sig_z)] <- 0

cat("\n== bulk (GSE51588) ==\n")
con <- gzfile("data/GSE51588_series_matrix.txt.gz", "rt"); raw <- readLines(con); close(con)
tstart <- grep("^!series_matrix_table_begin", raw); tend <- grep("^!series_matrix_table_end", raw)
titles <- strsplit(gsub("\"", "", gsub("^!Sample_title\t", "", raw[grep("^!Sample_title", raw)])), "\t")[[1]]
cat("n samples:", length(titles), "| examples:", paste(head(titles, 3), collapse=" | "), "\n")
body <- raw[(tstart + 2):(tend - 1)]; body <- body[!grepl("^!", body)]
hdr  <- strsplit(gsub("\"", "", raw[tstart + 1]), "\t")[[1]]
tab  <- read.delim(text = paste(body, collapse = "\n"), header = FALSE, quote = "",
                   stringsAsFactors = FALSE, check.names = FALSE)
colnames(tab) <- hdr
idcol <- grep("^ID_REF$|^ID_ref$", hdr, ignore.case = TRUE)[1]
cat("ID column:", hdr[idcol], "\n")
probe <- gsub("\"", "", as.character(tab[[idcol]]))
expr  <- as.matrix(tab[, setdiff(seq_along(hdr), idcol), drop = FALSE])
expr  <- matrix(suppressWarnings(as.numeric(gsub("\"", "", expr))), nrow = nrow(expr))
colnames(expr) <- hdr[setdiff(seq_along(hdr), idcol)]
cat("bulk matrix:", dim(expr), "\n")

p2s <- read.delim("data/GPL13497_probe2symbol.tsv", stringsAsFactors = FALSE)
sym <- p2s$symbol[match(probe, p2s$probe)]
keep <- !is.na(sym) & sym != ""; expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
cat("probes mapped:", nrow(expr), "\n")
tb <- table(sym)
expr_g <- rowsum(expr, group = sym, reorder = TRUE) / as.numeric(tb[rownames(rowsum(expr, group = sym, reorder = TRUE))])
cat("genes after collapsing:", nrow(expr_g), "\n")
if (max(expr_g, na.rm = TRUE) < 60) { cat("log2 -> linear\n"); expr_g <- 2^expr_g }

parts <- do.call(rbind, strsplit(titles, "-"))
disease <- parts[, 1]; comp <- parts[, 2]; donor <- parts[, 3]
cat("disease:", paste(unique(disease), collapse="/"), "| comp:", paste(unique(comp), collapse="/"),
    "| donors:", length(unique(donor)), "\n")

cat("\n== deconvolution ==\n")
common <- intersect(rownames(sig_z), rownames(expr_g)); cat("genes used:", length(common), "\n")
S <- sig_z[common, , drop = FALSE]
nnls_solve <- function(A, b) {
  r <- optim(rep(1/ncol(A), ncol(A)),
             function(x) sum((as.vector(A %*% x) - b)^2),
             function(x) as.vector(2 * t(A) %*% (as.vector(A %*% x) - b)),
             method = "L-BFGS-B", lower = rep(0, ncol(A)), control = list(maxit = 3000, factr = 1e7))
  x <- r$par; if (sum(x) > 0) x <- x / sum(x); x
}
props <- as.data.frame(t(sapply(seq_len(ncol(expr_g)), function(j) {
  vz <- as.vector(scale(expr_g[common, j])); vz[is.na(vz)] <- 0; nnls_solve(S, vz)
})))
colnames(props) <- subtypes
props$gsm <- colnames(expr_g); props$sample <- titles
props$disease <- disease; props$compartment <- comp; props$donor <- donor
cat("mean proportions by compartment:/n")
print(round(aggregate(props[, subtypes], by = list(comp = props$compartment), FUN = mean)[, -1], 4))
write.csv(props, file.path(RES, "33_deconvolution_proportions.csv"), row.names = FALSE)

paired_test <- function(df, s) {
  d <- intersect(df$donor[df$compartment == "LT"], df$donor[df$compartment == "MT"])
  a <- df[df$compartment == "LT" & df$donor %in% d, ]; a <- a[order(a$donor), s]
  b <- df[df$compartment == "MT" & df$donor %in% d, ]; b <- b[order(b$donor), s]
  tt <- t.test(a, b, paired = TRUE)
  data.frame(subtype = s, n_pairs = length(a), mean_LT = mean(a), mean_MT = mean(b),
             delta_LT_minus_MT = mean(a - b), t = tt$statistic, p = tt$p.value)
}
cat("\n== paired LT(non-load-bearing) vs MT(load-bearing): ALL ==\n")
paired <- do.call(rbind, lapply(subtypes, function(s) paired_test(props, s)))
paired$p_adj <- p.adjust(paired$p, "BH"); print(paired, digits = 4)
write.csv(paired, file.path(RES, "33_deconvolution_paired.csv"), row.names = FALSE)

cat("\n== paired, OA donors only ==\n")
oa <- props[props$disease == "OA", ]
pairedOA <- do.call(rbind, lapply(subtypes, function(s) paired_test(oa, s)))
pairedOA$p_adj <- p.adjust(pairedOA$p, "BH"); print(pairedOA, digits = 4)
write.csv(pairedOA, file.path(RES, "33_deconvolution_paired_OA.csv"), row.names = FALSE)

cat("\n== PIEZO1 vs proportions ==\n")
pz <- if ("PIEZO1" %in% rownames(expr_g)) expr_g["PIEZO1", ] else rep(NA_real_, ncol(expr_g))
props$PIEZO1 <- as.numeric(pz)
corr <- do.call(rbind, lapply(subtypes, function(s) {
  ct <- cor.test(props[[s]], props$PIEZO1, method = "spearman")
  data.frame(subtype = s, rho = ct$estimate, p = ct$p.value)
}))
corr$p_adj <- p.adjust(corr$p, "BH"); print(corr, digits = 4)
write.csv(corr, file.path(RES, "33_proportion_PIEZO1_corr.csv"), row.names = FALSE)

cat("\n== mechanotransduction (YAP/TAZ) ==\n")
yap <- c("CCN2","CTGF","CYR61","ANKRD1","AMOTL2","THBS1","NUAK2","LATS2","TEAD1","TEAD3",
         "AXIN2","CRIM1","IGFBP3","GADD45A","AMOTL1","F3","MYOF","PTX3")
yap <- intersect(yap, rownames(expr_g))
cat("genes present:", length(yap), ":", paste(yap, collapse=", "), "\n")
zexpr <- t(scale(t(log2(expr_g + 1))))
score <- colMeans(zexpr[yap, , drop = FALSE], na.rm = TRUE)
props$mech_score <- as.numeric(score)
d <- intersect(props$donor[props$compartment=="LT"], props$donor[props$compartment=="MT"])
a <- props[props$compartment=="LT" & props$donor %in% d, ]; a <- a[order(a$donor), ]
b <- props[props$compartment=="MT" & props$donor %in% d, ]; b <- b[order(b$donor), ]
tt <- t.test(a$mech_score, b$mech_score, paired = TRUE)
cat(sprintf("mech LT=%.4f MT=%.4f delta=%.4f t=%.3f p=%.4g\n",
            mean(a$mech_score), mean(b$mech_score), mean(a$mech_score-b$mech_score), tt$statistic, tt$p.value))
ct <- cor.test(props$mech_score, props$PIEZO1, method = "spearman")
cat(sprintf("mech vs PIEZO1: rho=%.4f p=%.4g\n", ct$estimate, ct$p.value))
write.csv(data.frame(comparison=c("LT_vs_MT_mech","mech_vs_PIEZO1"),
                     estimate=c(mean(a$mech_score-b$mech_score), ct$estimate),
                     p=c(tt$p.value, ct$p.value)),
          file.path(RES, "33_mechano_score.csv"), row.names = FALSE)
write.csv(props, file.path(RES, "33_sample_table.csv"), row.names = FALSE)
cat("\nDONE\n")
