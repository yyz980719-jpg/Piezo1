# Replication of the compartment gradient in single-cell cartilage (GSE152805)
# 3 OA donors x paired lateral (oLT, non-load-bearing) vs medial (MT, load-bearing) tibial cartilage
source("scripts/_bootstrap.R")
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
RES <- PIEZO1_RESULTS; set.seed(1)
dir <- data_path("GSE152805_raw")
files <- list.files(dir)
samples <- unique(sub("\\.(barcodes|genes|matrix).*", "", files[grep("barcodes", files, invert=FALSE)]))
# build sample list from matrix files
mtx <- files[grep("matrix.mtx.gz$", files)]
cat("samples found:", length(mtx), "\n")
keep <- c("OA_oLT_113","OA_oLT_116","OA_oLT_118","OA_MT_113","OA_MT_116","OA_MT_118")

counts <- NULL; meta <- data.frame()
for (s in keep) {
  gm <- sub("_OA.*","",s); gm <- sub("^GSM[0-9]+_","", s)
  f_m <- file.path(dir, list.files(dir, pattern=paste0("^GSM[0-9]+_", s, "\\.matrix")))
  f_b <- file.path(dir, list.files(dir, pattern=paste0("^GSM[0-9]+_", s, "\\.barcodes")))
  f_g <- file.path(dir, list.files(dir, pattern=paste0("^GSM[0-9]+_", s, "\\.genes")))
  M <- readMM(gzfile(f_m))
  bc <- read.delim(gzfile(f_b), header=FALSE, stringsAsFactors=FALSE)[,1]
  gn <- read.delim(gzfile(f_g), header=FALSE, stringsAsFactors=FALSE)
  dimnames(M) <- list(gn[,2], bc)
  cat(sprintf("  %s: %d genes x %d cells\n", s, nrow(M), ncol(M)))
  if (is.null(counts)) counts <- M else counts <- cbind(counts, M)
  comp <- if (grepl("oLT", s)) "lateral_nonload" else "medial_load"
  don  <- sub(".*_", "", s)
  meta <- rbind(meta, data.frame(cell=paste(s, bc, sep="_"), compartment=comp, donor=don))
}
rownames(meta) <- meta$cell
meta$compartment <- factor(meta$compartment, levels=c('medial_load','lateral_nonload'))
cat("combined:", dim(counts), "\n")

# QC
nc <- Matrix::colSums(counts); nf <- Matrix::colSums(counts > 0)
mt <- grep("^MT-", rownames(counts), value=TRUE)
mp <- if (length(mt)) Matrix::colSums(counts[mt,,drop=FALSE])/nc else rep(0, ncol(counts))
ok <- nf >= 200 & nc >= 500 & mp < 0.25
cat("cells passing QC:", sum(ok), "of", ncol(counts), "\n")
counts <- counts[, ok, drop=FALSE]; meta <- meta[ok, ]
# remove duplicated genes (keep first occurrence)
if (any(duplicated(rownames(counts)))) counts <- counts[!duplicated(rownames(counts)), ]

# log-normalize
lib <- Matrix::colSums(counts)
expr <- Matrix::t(log1p(Matrix::t(counts) / lib * 1e4))
expr <- as.matrix(expr)

# build programme signatures from GSE104782 (as in 33)
seu <- readRDS("data/scrna/GSE104782_seurat.rds")
sub <- as.character(seu@meta.data$subtype)
d0 <- tryCatch(GetAssayData(seu, layer="data"), error=function(e) GetAssayData(seu, slot="data"))
lin <- as.matrix(expm1(d0)); subtypes <- sort(unique(sub))
mean_mat <- sapply(subtypes, function(s) rowMeans(lin[, sub==s, drop=FALSE]))
sig_genes <- unique(unlist(lapply(subtypes, function(s){
  o <- rowMeans(mean_mat[, setdiff(subtypes,s), drop=FALSE])
  fc <- (mean_mat[,s]+1e-6)/(o+1e-6); fc[mean_mat[,s]<0.05] <- 0
  names(sort(fc, decreasing=TRUE))[1:100] })))
sig <- mean_mat[sig_genes,,drop=FALSE]; sig_z <- t(scale(t(sig))); sig_z[is.na(sig_z)] <- 0

# score GSE152805 cells
common <- intersect(sig_genes, rownames(expr))
zz <- t(scale(t(expr[common,,drop=FALSE]))); zz[is.na(zz)] <- 0
W <- t(zz) %*% sig_z[common,,drop=FALSE]          # cells x programmes
colnames(W) <- subtypes
pz <- if ("PIEZO1" %in% rownames(expr)) as.numeric(expr["PIEZO1",]) else rep(NA, ncol(expr))
cat("cells scored:", nrow(W), "| PIEZO1 detected in", round(mean(pz>0)*100,1), "% of cells\n")

# paired cell-level tests with donor blocking (linear model)
out <- do.call(rbind, lapply(subtypes, function(s){
  fit <- lm(W[,s] ~ meta$compartment + meta$donor)
  cf <- summary(fit)$coefficients
  data.frame(programme=s, delta_lateral_minus_medial=cf["meta$compartmentlateral_nonload","Estimate"],
             t=cf["meta$compartmentlateral_nonload","t value"], p=cf["meta$compartmentlateral_nonload","Pr(>|t|)"]) }))
out$p_adj <- p.adjust(out$p, "BH")
cat("\n== programme weights, lateral(non-load-bearing) vs medial(load-bearing) ==\n"); print(out, digits=4)

fitp <- lm(pz ~ meta$compartment + meta$donor)
cfp <- summary(fitp)$coefficients
cat(sprintf("\nPIEZO1 (single-cell, donor-blocked): delta(lateral-medial) = %.4f, P = %.3g\n",
            cfp["meta$compartmentlateral_nonload","Estimate"], cfp["meta$compartmentlateral_nonload","Pr(>|t|)"]))

# per-donor direction consistency (n = 3 donors -> direction is the honest statistic)
dons <- unique(meta$donor)
cat("\n== per-donor mean weight difference (lateral - medial) ==\n")
pdiff <- t(sapply(dons, function(dd){
  sapply(subtypes, function(s){
    a <- W[meta$donor==dd & meta$compartment=="lateral_nonload", s]
    b <- W[meta$donor==dd & meta$compartment=="medial_load", s]
    mean(a) - mean(b) }) }))
colnames(pdiff) <- subtypes; rownames(pdiff) <- dons
print(round(pdiff,1))
cat("donors with HomC higher in lateral:", sum(pdiff[,"HomC"]>0), "/ 3\n")
cat("donors with preHTC higher in medial:", sum(pdiff[,"preHTC"]<0), "/ 3\n")
cat("donors with EC higher in lateral:", sum(pdiff[,"EC"]>0), "/ 3\n")
cat("donors with ProC higher in medial:", sum(pdiff[,"ProC"]<0), "/ 3\n")

# pseudo-bulk PIEZO1 per donor
pbg <- t(sapply(dons, function(dd){
  sapply(c("lateral_nonload","medial_load"), function(cp){
    idx <- which(meta$donor==dd & meta$compartment==cp)
    mean(pz[idx]) }) }))
rownames(pbg) <- dons; colnames(pbg) <- c("lateral","medial")
cat("\nmean PIEZO1 log-CPM per cell, lateral vs medial, per donor:\n"); print(round(pbg,3))

write.csv(out, file.path(RES,"39_gse152805_programmes.csv"), row.names=FALSE)
# Keep programme differences and PIEZO1 donor summaries in separate rectangular
# tables.  The historical output appended a three-donor PIEZO1 vector to a
# seven-programme matrix, which silently recycled values across columns.
write.csv(pdiff, file.path(RES,"39_gse152805_perdonor.csv"))
write.csv(pbg, file.path(RES,"39_gse152805_piezo1.csv"))
cat("\nDONE\n")
