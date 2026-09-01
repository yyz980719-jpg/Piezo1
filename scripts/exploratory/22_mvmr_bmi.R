#!/usr/bin/env Rscript
# 22_mvmr_bmi.R -- Phase 2-B: multivariable MR to strip BMI confounding
# from the TF -> OA effect (B design). Uses:
#   - exposure 1 : TF cis-eQTL (beta_x, se_x)         [from 17 multi_snp]
#   - confounder : BMI (finngen_R12_BMI_IRN.gz)         [pulled from GCS]
#   - outcome    : OA  (beta_y, se_y)                  [from 17 multi_snp]
# All three are oriented to the FinnGen OA effect allele (fg_alt), so they
# are mutually allele-harmonized. R is launched via /d/R/rscript.sh to avoid
# the Git-Bash POSIX-path segfault.
suppressMessages({ library(data.table) })

source("scripts/_bootstrap.R")
ROOT  <- PIEZO1_ROOT
ms    <- fread(file.path(ROOT,"results","17_MR_TF_OA_multi_snp.csv"),
               colClasses=list(character=c("SNP","tf","outcome","fg_alt","fg_ref","tissue")))
bmi   <- fread(file.path(ROOT,"results","_bmi_subset.tsv"))   # has BMI header
setnames(bmi, names(bmi), c("chrom","pos","ref","alt","rsids","gene","pval","mlogp","beta","sebeta","af_alt"))

# exact-token match of our SNPs in the BMI rsids field
snps <- unique(ms$SNP)
bmi[, rs_tokens := strsplit(rsids, " ")]
bmi[, keep := FALSE]
for (s in snps) bmi[ sapply(rs_tokens, function(v) s %in% v), keep := TRUE ]
bmi <- bmi[ keep==TRUE ]
# collapse to one row per SNP (keep first match; multi-matches are rare, same variant)
bmi <- bmi[ , .SD[1], by="rsids" ]   # not ideal key; instead key by token below

# robust: build SNP->BMI lookup via token scan
bmi_lu <- list()
for (i in seq_len(nrow(bmi))) {
  toks <- bmi$rs_tokens[[i]]
  for (t in toks) if (t %in% snps) bmi_lu[[ t ]] <- list(alt=bmi$alt[i], ref=bmi$ref[i],
                                                          beta=bmi$beta[i], se=bmi$sebeta[i])
}
cat(sprintf("BMI lookup built for %d SNPs\n", length(bmi_lu)))

# orient BMI beta to fg_alt (FinnGen OA effect allele)
ms[, b_bmi := NA_real_]
ms[, se_bmi := NA_real_]
for (i in seq_len(nrow(ms))) {
  s <- ms$SNP[i]; L <- bmi_lu[[s]]; if (is.null(L)) next
  if (L$alt == ms$fg_alt[i])      { ms$b_bmi[i] <-  L$beta; ms$se_bmi[i] <- L$se }
  else if (L$alt == ms$fg_ref[i]) { ms$b_bmi[i] <- -L$beta; ms$se_bmi[i] <- L$se }
  # else: strand/asm mismatch -> leave NA
}
ms2 <- ms[ !is.na(b_bmi) & !is.na(beta_x) & !is.na(beta_y) ]
cat(sprintf("rows with TF + BMI + OA harmonized: %d (of %d)\n", nrow(ms2), nrow(ms)))

# ---- MVMR (Sanderson & Davey Smith 2019): GLS with weights 1/se_y^2 ----
mvmr <- function(bx, sbx, bb, sbb, by, sby) {
  n <- length(by)
  if (n < 3) return(list(n=n, ok=FALSE))
  X <- cbind(bx, bb)                       # n x 2 exposures
  Y <- by
  W <- diag(1 / (sby^2))
  XtW  <- t(X) %*% W
  XtWX <- XtW %*% X
  inv  <- tryCatch(solve(XtWX), error=function(e) NULL)
  if (is.null(inv)) return(list(n=n, ok=FALSE, singular=TRUE))
  G    <- inv %*% XtW %*% Y               # 2x1
  Cov  <- inv
  se   <- sqrt(diag(Cov))
  # unadjusted IVW (TF only)
  x  <- matrix(bx, ncol=1)
  i0 <- solve(t(x) %*% W %*% x)
  g0 <- as.numeric(i0 %*% t(x) %*% W %*% Y)
  se0<- sqrt(as.numeric(i0))
  # Q (residual heterogeneity)
  res <- Y - X %*% G
  Q   <- as.numeric(t(res) %*% W %*% res); Qdf <- n - 2
  Qp  <- pchisq(Q, Qdf, lower.tail=FALSE)
  # modified Q per exposure (drop column j)
  mq <- rep(NA_real_, 2)
  for (j in 1:2) {
    Xd <- X[, -j, drop=FALSE]
    invd <- tryCatch(solve(t(Xd)%*%W%*%Xd), error=function(e) NULL)
    if (is.null(invd)) next
    Gd <- invd %*% t(Xd) %*% W %*% Y
    rd <- Y - Xd %*% Gd
    Qd <- as.numeric(t(rd) %*% W %*% rd)
    mq[j] <- Qd - Q
  }
  mqp <- sapply(mq, function(q) if(is.na(q)) NA else pchisq(q,1,lower.tail=FALSE))
  list(n=n, ok=TRUE,
       b_TF_adj=G[1], se_TF_adj=se[1], p_TF_adj=2*pnorm(-abs(G[1]/se[1])),
       b_BMI_adj=G[2], se_BMI_adj=se[2], p_BMI_adj=2*pnorm(-abs(G[2]/se[2])),
       b_TF_unadj=g0, se_TF_unadj=se0, p_TF_unadj=2*pnorm(-abs(g0/se0)),
       Q=Q, Qdf=Qdf, Qp=Qp, modQ_TF=mq[1], modQ_TF_p=mqp[1],
       modQ_BMI=mq[2], modQ_BMI_p=mqp[2])
}
out <- list()
for (oc_i in unique(ms2$outcome)) for (tf_i in unique(ms2$tf)) {
  sub <- ms2[ tf==tf_i & outcome==oc_i ]
  if (nrow(sub)==0) next
  r <- mvmr(sub$beta_x, sub$se_x, sub$b_bmi, sub$se_bmi, sub$beta_y, sub$se_y)
  out[[ sprintf("%s|%s", tf_i, oc_i) ]] <- data.table(
    tf=tf_i, outcome=oc_i, n_snp=r$n,
    TF_unadj_beta=if(r$ok) round(r$b_TF_unadj,4) else NA,
    TF_unadj_se =if(r$ok) round(r$se_TF_unadj,4) else NA,
    TF_unadj_p  =if(r$ok) signif(r$p_TF_unadj,3) else NA,
    TF_adj_beta =if(r$ok) round(r$b_TF_adj,4) else NA,
    TF_adj_se   =if(r$ok) round(r$se_TF_adj,4) else NA,
    TF_adj_p    =if(r$ok) signif(r$p_TF_adj,3) else NA,
    BMI_adj_beta=if(r$ok) round(r$b_BMI_adj,4) else NA,
    BMI_adj_se  =if(r$ok) round(r$se_BMI_adj,4) else NA,
    BMI_adj_p   =if(r$ok) signif(r$p_BMI_adj,3) else NA,
    Q_total     =if(r$ok) round(r$Q,3) else NA,
    Q_df        =if(r$ok) r$Qdf else NA,
    Q_p         =if(r$ok) signif(r$Qp,3) else NA,
    modQ_TF_p   =if(r$ok) signif(r$modQ_TF_p,3) else NA,
    modQ_BMI_p  =if(r$ok) signif(r$modQ_BMI_p,3) else NA,
    mvmr_ok     =r$ok)
}
res <- rbindlist(out, fill=TRUE)
fwrite(res, file.path(ROOT,"results","22_mvmr_bmi.csv"))
cat("=== MVMR-BMI results ===\n"); print(res, nrow=Inf)

# ---- figure: TF->OA effect unadjusted vs BMI-adjusted ----
library(ggplot2)
pl <- res[ !is.na(TF_adj_beta) ]
pl[, lab := paste0(tf,"/",outcome)]
pl[, ymin_u := TF_unadj_beta-1.96*TF_unadj_se]
pl[, ymax_u := TF_unadj_beta+1.96*TF_unadj_se]
pl[, ymin_a := TF_adj_beta-1.96*TF_adj_se]
pl[, ymax_a := TF_adj_beta+1.96*TF_adj_se]
png(file.path(ROOT,"figs","fig71_mvmr_bmi.png"), width=900, height=600, res=110)
# forest: two rows per tf/outcome (unadj, adj)
ord <- order(pl$TF_unadj_beta)
d <- pl[ord]
par(mfrow=c(1,1), mar=c(4,9,3,1))
y <- seq_len(nrow(d))
plot(d$TF_unadj_beta, y, xlim=range(c(d$ymin_u,d$ymax_u,d$ymin_a,d$ymax_a)),
     yaxt="n", ylab="", xlab="TF -> OA effect (IVW beta, allele-harmonized)", pch=16, col=2)
arrows(d$ymin_u, y, d$ymax_u, y, len=0.05, col=2, code=3)
points(d$TF_adj_beta, y+0.18, pch=17, col=4)
arrows(d$ymin_a, y+0.18, d$ymax_a, y+0.18, len=0.05, col=4, code=3)
abline(v=0, lty=2)
axis(2, at=y, labels=d$lab, las=2, cex=0.7)
legend("topright", c("unadjusted","BMI-adjusted"), pch=c(16,17), col=c(2,4))
title("MVMR-BMI: TF->OA effect before vs after BMI adjustment")
dev.off()
cat("fig71 written\n")
