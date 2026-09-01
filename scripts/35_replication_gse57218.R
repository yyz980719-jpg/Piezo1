# Replication of the loading/programme gradient in an independent, paired CARTILAGE dataset (GSE57218)
source("scripts/_bootstrap.R")
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
RES <- PIEZO1_RESULTS; set.seed(1)

cat("== GSE57218 (paired preserved vs OA-affected cartilage, RAAK) ==\n")
con <- gzfile("data/GSE57218_series_matrix.txt.gz","rt"); raw <- readLines(con); close(con)
tstart <- grep("^!series_matrix_table_begin", raw); tend <- grep("^!series_matrix_table_end", raw)
titles <- strsplit(gsub("\"","",gsub("^!Sample_title\t","",raw[grep("^!Sample_title",raw)])),"\t")[[1]]
body <- raw[(tstart+2):(tend-1)]; body <- body[!grepl("^!",body)]
hdr <- strsplit(gsub("\"","",raw[tstart+1]),"\t")[[1]]
tab <- read.delim(text=paste(body,collapse="\n"), header=FALSE, quote="", stringsAsFactors=FALSE, check.names=FALSE)
colnames(tab) <- hdr
idcol <- grep("^ID_REF$|^ID$", hdr, ignore.case=TRUE)[1]
probe <- gsub("\"","",as.character(tab[[idcol]]))
expr <- as.matrix(tab[, setdiff(seq_along(hdr), idcol), drop=FALSE])
expr <- matrix(suppressWarnings(as.numeric(gsub("\"","",expr))), nrow=nrow(expr))
colnames(expr) <- hdr[setdiff(seq_along(hdr), idcol)]
cat("matrix:", dim(expr), "\n")

# probe -> symbol (GPL6947 annot)
zf <- gzfile("data/GPL6947.annot.gz","rt"); an <- readLines(zf); close(zf)
an <- an[!grepl("^\\^|^!", an)]
ah <- which(grepl("^ID\t", an))[1]
if (is.na(ah)) ah <- which(grepl("^ID_REF\t|^Probe", an))[1]
ahdr <- strsplit(an[ah],"\t")[[1]]
i_id <- grep("^ID$|^ID_REF$|^Probe", ahdr)[1]
i_sym<- grep("^Gene symbol$|^Symbol$", ahdr)[1]
cat("annot cols:", ahdr[i_id], "/", ahdr[i_sym], "\n")
sp <- strsplit(an[(ah+1):length(an)], "\t", fixed=TRUE)
pid <- sapply(sp, function(x) if (length(x) >= max(i_id,i_sym)) trimws(x[i_id]) else NA)
psy <- sapply(sp, function(x) if (length(x) >= max(i_id,i_sym)) trimws(x[i_sym]) else NA)
p2s <- data.frame(probe=pid, symbol=psy, stringsAsFactors=FALSE)
p2s <- p2s[!is.na(p2s$symbol) & p2s$symbol != "" & p2s$symbol != "NA", ]
sym <- p2s$symbol[match(probe, p2s$probe)]
keep <- !is.na(sym) & sym!=""
expr <- expr[keep,,drop=FALSE]; sym <- sym[keep]
tb <- table(sym)
expr_g <- rowsum(expr, group=sym, reorder=TRUE) / as.numeric(tb[rownames(rowsum(expr, group=sym, reorder=TRUE))])
cat("genes:", nrow(expr_g), "\n")

# metadata
raak <- sub("cartilage_RAAK_", "", titles); raak <- sub("_.*$", "", raak)
ttype <- sub(".*_", "", titles)
cat("types:", paste(names(table(ttype)), table(ttype), sep="=", collapse=" "), "\n")

# ---------- 1. PIEZO1 paired test ----------
cat("\n== PIEZO1: preserved vs OA (paired) ==\n")
pz <- if ("PIEZO1" %in% rownames(expr_g)) expr_g["PIEZO1",] else rep(NA_real_, ncol(expr_g))
res <- data.frame()
pa <- function(gene){
  v <- if (gene %in% rownames(expr_g)) expr_g[gene,] else rep(NA_real_, ncol(expr_g))
  ids <- intersect(raak[ttype=="Preserved"], raak[ttype=="OA"])
  if (all(is.na(v))) return(data.frame(gene=gene, n_pairs=0, mean_preserved=NA, mean_OA=NA,
      delta_preserved_minus_OA=NA, t=NA, p=NA, OA_vs_healthy_p=NA))
  a <- v[ttype=="Preserved" & raak %in% ids]; b <- v[ttype=="OA" & raak %in% ids]
  na <- raak[ttype=="Preserved" & raak %in% ids]; nb <- raak[ttype=="OA" & raak %in% ids]
  a <- a[order(na)]; b <- b[order(nb)]
  tt <- t.test(a, b, paired=TRUE)
  h <- v[ttype=="Healthy"]
  tt2 <- if (length(h)>1) t.test(v[ttype=="OA"], h) else NULL
  data.frame(gene=gene, n_pairs=length(a), mean_preserved=mean(a), mean_OA=mean(b),
             delta_preserved_minus_OA=mean(a-b), t=tt$statistic, p=tt$p.value,
             OA_vs_healthy_p=ifelse(is.null(tt2), NA, tt2$p.value))
}
res <- rbind(res, pa("PIEZO1"))
for (g in c("COL2A1","SOX9","ACAN","PRG4","MMP13","RUNX2","IL1B","TNF")) res <- rbind(res, pa(g))
print(res, digits=4)

# ---------- 2. chondrocyte programme weights (GSE104782 signatures) ----------
cat("\n== chondrocyte programme weights: preserved vs OA ==\n")
seu <- readRDS("data/scrna/GSE104782_seurat.rds"); sub <- as.character(seu@meta.data$subtype)
dat <- tryCatch(GetAssayData(seu, layer="data"), error=function(e) GetAssayData(seu, slot="data"))
lin <- as.matrix(expm1(dat)); subtypes <- sort(unique(sub))
mean_mat <- sapply(subtypes, function(s) rowMeans(lin[, sub==s, drop=FALSE]))
sig_genes <- unique(unlist(lapply(subtypes, function(s){
  others <- rowMeans(mean_mat[, setdiff(subtypes,s), drop=FALSE])
  fc <- (mean_mat[,s]+1e-6)/(others+1e-6); fc[mean_mat[,s]<0.05] <- 0
  names(sort(fc, decreasing=TRUE))[1:100]
})))
sig <- mean_mat[sig_genes,,drop=FALSE]; sig_z <- t(scale(t(sig))); sig_z[is.na(sig_z)] <- 0
common <- intersect(rownames(sig_z), rownames(expr_g)); S <- sig_z[common,,drop=FALSE]
nnls_solve <- function(A,b){
  r <- optim(rep(1/ncol(A),ncol(A)), function(x) sum((as.vector(A%*%x)-b)^2),
             function(x) as.vector(2*t(A)%*%(as.vector(A%*%x)-b)),
             method="L-BFGS-B", lower=rep(0,ncol(A)), control=list(maxit=3000,factr=1e7))
  x <- r$par; if (sum(x)>0) x <- x/sum(x); x }
W <- t(sapply(seq_len(ncol(expr_g)), function(j){
  vz <- as.vector(scale(expr_g[common,j])); vz[is.na(vz)] <- 0; nnls_solve(S, vz) }))
colnames(W) <- subtypes
ids <- intersect(raak[ttype=="Preserved"], raak[ttype=="OA"])
pw <- do.call(rbind, lapply(subtypes, function(s){
  a <- W[ttype=="Preserved" & raak %in% ids, s]; b <- W[ttype=="OA" & raak %in% ids, s]
  na <- raak[ttype=="Preserved" & raak %in% ids]; nb <- raak[ttype=="OA" & raak %in% ids]
  a <- a[order(na)]; b <- b[order(nb)]; tt <- t.test(a,b,paired=TRUE)
  data.frame(subtype=s, mean_preserved=mean(a), mean_OA=mean(b),
             delta=mean(a-b), t=tt$statistic, p=tt$p.value) }))
pw$p_adj <- p.adjust(pw$p,"BH"); print(pw, digits=4)

# ---------- 3. transfer compartment signatures from GSE51588 ----------
cat("\n== transferred compartment signatures (GSE51588 -> GSE57218) ==\n")
degf <- file.path(RES,"07_DEG_paired_MedialVsLateral_OA.csv")
if (file.exists(degf)) {
  d <- read.csv(degf, stringsAsFactors=FALSE)
  cat("DEG cols:", paste(head(names(d),12), collapse=", "), "\n")
  gcol <- names(d)[grep("^gene|^Gene|^SYMBOL|^symbol", names(d))][1]
  lcol <- names(d)[grep("logFC|log2FC|log2fc", names(d))][1]
  if (!is.na(gcol) && !is.na(lcol)) {
    d <- d[!is.na(d[[lcol]]), ]
    # logFC is Medial - Lateral, so POSITIVE logFC = higher in MEDIAL (load-bearing),
    #                            NEGATIVE logFC = higher in LATERAL (non-load-bearing)
    pz_fc <- d[[lcol]][match("PIEZO1", d[[gcol]])]
    cat("PIEZO1 logFC (Medial-Lateral) =", pz_fc, " -> higher in",
        ifelse(!is.na(pz_fc) && pz_fc < 0, "LATERAL (non-load-bearing)", "MEDIAL"), "\n")
    up_med <- head(d[order(-d[[lcol]]), gcol], 200)   # load-bearing (medial) genes
    up_lat <- head(d[order(d[[lcol]]), gcol], 200)    # non-load-bearing (lateral) genes
    cat("lateral-up(non-loaded) genes:", length(up_lat), " medial-up(loaded):", length(up_med), "\n")
    zx <- t(scale(t(expr_g)))
    score <- function(gs){
      gs <- intersect(gs, rownames(zx))
      if (length(gs) < 5) return(rep(NA_real_, ncol(zx)))
      colMeans(zx[gs,,drop=FALSE], na.rm=TRUE) }
    s_lat <- score(up_lat); s_med <- score(up_med)
    a1 <- s_lat[ttype=="Preserved" & raak %in% ids]; b1 <- s_lat[ttype=="OA" & raak %in% ids]
    n1 <- raak[ttype=="Preserved" & raak %in% ids]; n2 <- raak[ttype=="OA" & raak %in% ids]
    a1 <- a1[order(n1)]; b1 <- b1[order(n2)]
    t1 <- t.test(a1,b1,paired=TRUE)
    a2 <- s_med[ttype=="Preserved" & raak %in% ids]; b2 <- s_med[ttype=="OA" & raak %in% ids]
    a2 <- a2[order(n1)]; b2 <- b2[order(n2)]
    t2 <- t.test(a2,b2,paired=TRUE)
    cat(sprintf("non-load-bearing(lateral-up) signature: preserved-OA delta=%.4f t=%.3f p=%.4g\n",
                mean(a1-b1), t1$statistic, t1$p.value))
    cat(sprintf("load-bearing(medial-up) signature:     preserved-OA delta=%.4f t=%.3f p=%.4g\n",
                mean(a2-b2), t2$statistic, t2$p.value))
    write.csv(data.frame(signature=c("lateral_up_nonloaded","medial_up_loaded"),
                         delta_preserved_minus_OA=c(mean(a1-b1), mean(a2-b2)),
                         t=c(t1$statistic,t2$statistic), p=c(t1$p.value,t2$p.value)),
              file.path(RES,"35_replication_signatures.csv"), row.names=FALSE)
  } else cat("DEG gene/logFC columns not found\n")
} else cat("DEG file not found\n")

write.csv(res, file.path(RES,"35_replication_genes.csv"), row.names=FALSE)
write.csv(pw,  file.path(RES,"35_replication_programmes.csv"), row.names=FALSE)
cat("\nDONE\n")
