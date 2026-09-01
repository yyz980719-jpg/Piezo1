# =====================================================================
# 17_mr_tf_oa.R  MR: 8 consensus TF expression -> OA risk (FinnGen R12)
# Exposure: GTEx tissue eQTLs (Whole Blood + connective tissues),
#           significant cis-eQTLs as instruments (p<5e-8, fallback p<1e-5)
# Outcome: FinnGen R12 M13_ARTHROSIS / _KNEE / _COX (gz in C: workspace)
# MR: Wald (1-SNP), IVW(RE), weighted median, MR-Egger, leave-one-out, Steiger
# Output: results/17_MR_TF_OA_*.csv ; figs/fig64_*.png
# =====================================================================
source("scripts/_bootstrap.R")
setwd(data_path("mr"))
suppressPackageStartupMessages({ library(gtexr); library(data.table); library(ggplot2) })
options(gtexr.itemsPerPage = 100000)
options(error = function() { cat("\n*** ERROR TRAP ***\n"); traceback(2); quit(status=1, save="no") })
dir_res <- PIEZO1_RESULTS
dir_fig <- PIEZO1_FIGURES

tfs <- c("TCF7L1","GATAD2A","PRDM16","BCL6","ZNF92","ZNF853","PKNOX2","ATF6")
tissue_order <- c("Whole_Blood","Muscle_Skeletal","Skin_Not_Sun_Exposed_Suprapubic",
                 "Artery_Aorta","Adipose_Subcutaneous")
finngen_dir <- data_path("mr")
gz <- c(ART = file.path(finngen_dir,"finngen_R12_M13_ARTHROSIS.gz"),
         KNEE = file.path(finngen_dir,"finngen_R12_M13_ARTHROSIS_KNEE.gz"),
         COX  = file.path(finngen_dir,"finngen_R12_M13_ARTHTROSIS_COX.gz"))
pheno_n <- list(ART=c(101454,315115), KNEE=c(61356,315115), COX=c(30802,315115))
plink  <- data_path("mr", if (.Platform$OS.type == "windows") "plink.exe" else "plink")
ref    <- data_path("mr", "ref_eur")

## ---------- MR core (reused from 11_mr_piezo1.R) ----------
ivw_re <- function(r, se) {
  w <- 1/se^2; k <- length(r); if (k==1) return(c(est=r, se=se))
  f <- sum(w*r)/sum(w); Q <- sum(w*(r-f)^2)
  tau2 <- max(0, (Q-k+1)/(sum(w)-sum(w^2)/sum(w)))
  wre <- 1/(se^2+tau2); est <- sum(wre*r)/sum(wre); s <- sqrt(1/sum(wre))
  c(est=est, se=s, Q=Q, df=k-1, tau2=tau2)
}
wmedian <- function(r, se, B=5000) {
  w <- 1/se^2; o <- order(r); cw <- cumsum(w[o])/sum(w)
  est <- r[o][which(cw>=0.5)[1]]
  rb <- replicate(B, { rj <- rnorm(length(r), r, se)
    oo <- order(rj); cww <- cumsum(w[oo])/sum(w); rj[oo][which(cww>=0.5)[1]] })
  c(est=est, se=sd(rb, na.rm=TRUE))
}
mr_egger <- function(bx, by, se_y) {
  fit <- summary(lm(by ~ bx, weights=1/se_y^2))
  c(slope=coef(fit)[2,1], se_slope=coef(fit)[2,2], intercept=coef(fit)[1,1],
    se_int=coef(fit)[1,2], p_int=coef(fit)[1,4])
}

## ---------- 1. resolve gencodeIds ----------
resolve <- function(sym){
  gs <- tryCatch(get_gene_search(geneId=sym), error=function(e) NULL)
  if(is.null(gs)||nrow(gs)==0) return(NA_character_)
  hit <- gs[tolower(gs$geneSymbol)==tolower(sym) & !is.na(gs$entrezGeneId),]
  if(nrow(hit)==0) hit <- gs[tolower(gs$geneSymbol)==tolower(sym),]
  if(nrow(hit)==0) return(NA_character_)
  as.character(hit$gencodeId[1])
}
cat("=== resolving TF gencodeIds ===\n")
gid <- sapply(tfs, resolve); print(gid)

## ---------- 2. fetch GTEx eQTLs (cached) ----------
cache <- data_path("mr", "17_gtex_eqtl_cache.rds")
if(file.exists(cache)){
  eqall <- readRDS(cache); cat("loaded cached GTEx eQTLs\n")
} else {
  eqall <- list()
  for(s in tfs){
    g <- gid[s]; if(is.na(g)){ eqall[[s]]<-NULL; next }
    best_tis <- NA; best_mp <- Inf; best_eq <- NULL
    per_tis <- list()
    for(tis in tissue_order){
      eq0 <- tryCatch(get_significant_single_tissue_eqtls(gencodeIds=g, tissueSiteDetailIds=tis),
                      error=function(e) NULL)
      if(is.null(eq0)||nrow(eq0)==0) next
      eq <- as.data.table(eq0)
      eq$pv <- as.numeric(as.character(eq$pValue))
      eq$nesv <- as.numeric(as.character(eq$nes))
      eq$tissue <- tis
      per_tis[[tis]] <- eq
      mp <- min(eq$pv, na.rm=TRUE)
      if(mp < best_mp){ best_mp<-mp; best_tis<-tis; best_eq<-eq }
    }
    eqall[[s]] <- list(gencodeId=g, tissue=best_tis, min_p=best_mp,
                       rows=if(is.null(best_eq)) NULL else best_eq, all_tissues=per_tis)
    cat(sprintf("  %-8s tissue=%-35s min_p=%.2e n_sig=%d\n", s, best_tis, best_mp,
                ifelse(is.null(best_eq),0,nrow(best_eq))))
  }
  saveRDS(eqall, cache)
}

## ---------- 3. build instruments per TF ----------
## Primary instrument = the single most-significant cis-eQTL SNP per TF (8 independent IVs).
## Sensitivity instrument = all significant cis-eQTLs per TF (LD clumping NOT possible:
##   the local 1000G EUR reference panel lacks coverage of these eQTL SNPs).
instr <- list(); lead_instr <- list()
for(s in tfs){
  e <- eqall[[s]]
  if(is.null(e) || is.null(e$rows) || nrow(e$rows)==0){ next }
  d <- as.data.table(e$rows)
  d <- d[!is.na(pv) & !is.na(nesv)]
  sel <- d[pv < 5e-8]; thr <- 5e-8          # instrument threshold
  if(nrow(sel)==0){ sel <- d[pv < 1e-5]; thr <- 1e-5 }
  if(nrow(sel)==0) next
  vt <- strsplit(sel$variantId, "_")         # chr, pos, ref, alt
  ref_g <- sapply(vt, `[`, 3); alt_g <- sapply(vt, `[`, 4); pos <- as.integer(sapply(vt, `[`, 2))
  sel[, `:=`(rsid=snpId, ref_g=ref_g, alt_g=alt_g, pos=pos,
             beta_x=nesv, se_x=abs(nesv)/(-qnorm(pv/2)),
             p_x=pv, tissue=e$tissue, tf=s, thr=thr)]
  instr[[s]] <- sel[, .(tf, tissue, rsid, ref_g, alt_g, pos, beta_x, se_x, p_x, variantId, nesv, pv, thr)]
  l <- sel[which.min(p_x)]                   # lead (most significant) cis-eQTL
  lead_instr[[s]] <- l[, .(tf, tissue, rsid, ref_g, alt_g, pos, beta_x, se_x, p_x, variantId, nesv, pv, thr)]
  cat(sprintf("  %-8s n_sig=%d, lead=%s p=%.2e (tissue=%s)\n", s, nrow(sel), l$rsid, l$p_x, e$tissue))
}
if(length(instr)==0){ cat("NO TF has instruments -> cannot run MR.\n"); quit(save="no") }

## ---------- 4. LD clump: skipped (reference panel unavailable) ----------
## The local ref_eur panel (v_<pos>_<ref> IDs, ~235k SNPs) does not contain the GTEx
## eQTL SNPs, so PLINK --clump cannot match them. We therefore use the lead-SNP-per-TF
## design as the primary (independent by construction) and report all-cis-eQTL as sensitivity.

## ---------- 5. harmonize + MR (explode rsids, keyed join) ----------
run_mr <- function(instr_set){
  all_snp <- list(); main_res <- list(); loo_res <- list(); steer_res <- list()
  for(oc in names(gz)){
    fg <- fread(gz[oc], showProgress=FALSE)
    setnames(fg, make.names(names(fg)))   # #chrom -> X.chrom
    ncase <- pheno_n[[oc]][1]; nctrl <- pheno_n[[oc]][2]
    target <- unique(unlist(lapply(instr_set, function(x) x$rsid)))
    # fast candidate retrieval: first/second rsid per row, keyed binary search
    # (avoids a 21M-row grepl, >5 min per phenotype)
    fg[, rs1 := sub("([^;,: ]+).*", "\\1", as.character(rsids))]
    fg[, rs2 := sub("^[^;,: ]+[;,: ]+([^;,: ]+).*", "\\1", as.character(rsids))]
    setkey(fg, rs1, rs2)
    cand <- fg[rs1 %in% target | rs2 %in% target]
    tg <- cand[, .(rsid = unlist(strsplit(as.character(rsids), "[;, ]+"))),
              by = .(ref, alt, beta, sebeta, pval, af_alt)]
    setkey(tg, rsid)
    cat(sprintf("  [%s] target rsids=%d, matched in FinnGen=%d\n", oc, length(target),
                sum(unique(tg$rsid) %in% target)))
    for(s in names(instr_set)){
      d <- copy(instr_set[[s]])
      rows <- list()
      for(i in 1:nrow(d)){
        ri <- tg[rsid == d$rsid[i]]            # fast exact join
        if(nrow(ri)==0) next
        ri <- ri[1]
        fg_alt <- ri$alt; fg_ref <- ri$ref
        ag <- d$alt_g[i]; rg <- d$ref_g[i]
        if(!(ag %in% c(fg_alt,fg_ref)) || !(rg %in% c(fg_alt,fg_ref))) next   # strand/asm mismatch
        if(setequal(c(ag,rg), c("A","T")) || setequal(c(ag,rg), c("C","G"))) next  # palindromic
        bx <- d$beta_x[i]; sx <- d$se_x[i]
        if(ag == fg_alt)      bx_a <- bx
        else if(ag == fg_ref) bx_a <- -bx
        else next
        by <- ri$beta; sy <- ri$sebeta; py <- ri$pval; af <- ri$af_alt
        wald <- by/bx_a; se_w <- sqrt(sy^2/bx_a^2 + (by*sx)^2/bx_a^4)
        rows[[length(rows)+1]] <- data.table(tf=s, outcome=oc, tissue=d$tissue[i], SNP=d$rsid[i],
          fg_alt=fg_alt, fg_ref=fg_ref, alt_g=ag, beta_x=bx_a, se_x=sx,
          beta_y=by, se_y=sy, p_y=py, eaf=af, wald=wald, se_wald=se_w)
      }
      if(length(rows)==0) next
      mm <- rbindlist(rows)
      k <- nrow(mm)
      iv <- ivw_re(mm$wald, mm$se_wald)
      wm <- if(k>=3) wmedian(mm$wald, mm$se_wald) else c(est=NA,se=NA)
      eg <- if(k>=3) mr_egger(mm$beta_x, mm$beta_y, mm$se_y) else rep(NA_real_,5)
      main_res[[paste(oc,s)]] <- data.table(
        tf=s, outcome=oc, tissue=mm$tissue[1], n_snp=k, thr=d$thr[1],
        method=c("IVW (RE)","Weighted median","MR-Egger slope","MR-Egger intercept"),
        estimate=c(iv["est"], wm["est"], eg["slope"], eg["intercept"]),
        se=c(iv["se"], wm["se"], eg["se_slope"], eg["se_int"]),
        p=2*pnorm(abs(c(iv["est"],wm["est"],eg["slope"],eg["intercept"]))/
                  c(iv["se"],wm["se"],eg["se_slope"],eg["se_int"]), lower.tail=FALSE),
        Q=c(iv["Q"],NA,NA,NA), p_het=c(pchisq(iv["Q"],iv["df"],lower.tail=FALSE),NA,NA,NA))
      if(k>=2) for(i in 1:k){
        idx <- setdiff(1:k,i)
        lo <- if(length(idx)>=1) ivw_re(mm$wald[idx], mm$se_wald[idx])["est"] else NA
        loo_res[[length(loo_res)+1]] <- data.table(tf=s, outcome=oc, excluded=mm$SNP[i], est=unname(lo))
      }
      Zy <- mm$beta_y/mm$se_y; Ny <- ncase+nctrl
      r2y <- Zy^2/(Zy^2+Ny); r2x <- mm$beta_x^2/(mm$beta_x^2 + (1/mm$se_x^2))  # approx r2 from se
      steer_res[[paste(oc,s)]] <- data.table(tf=s, outcome=oc, SNP=mm$SNP, r2_exposure=r2x, r2_outcome=r2y,
                                            correct_direction = r2x>r2y)
      all_snp[[paste(oc,s)]] <- mm
    }
    rm(fg, cand, tg); gc()
  }
  list(snp=rbindlist(all_snp, fill=TRUE), main=rbindlist(main_res, fill=TRUE),
       loo=rbindlist(loo_res, fill=TRUE), steer=rbindlist(steer_res, fill=TRUE))
}

cat("\n=== PRIMARY: lead cis-eQTL SNP per TF (8 independent IVs) ===\n")
lead_out  <- run_mr(lead_instr)
cat("\n=== SENSITIVITY: all significant cis-eQTLs (LD-clump unavailable) ===\n")
multi_out <- run_mr(instr)

fwrite(lead_out$snp,   file.path(dir_res,"17_MR_TF_OA_lead_snp.csv"))
fwrite(lead_out$main,  file.path(dir_res,"17_MR_TF_OA_lead_main.csv"))
fwrite(lead_out$loo,   file.path(dir_res,"17_MR_TF_OA_lead_leaveoneout.csv"))
fwrite(lead_out$steer, file.path(dir_res,"17_MR_TF_OA_lead_steiger.csv"))
fwrite(multi_out$snp,  file.path(dir_res,"17_MR_TF_OA_multi_snp.csv"))
fwrite(multi_out$main, file.path(dir_res,"17_MR_TF_OA_multi_main.csv"))
fwrite(multi_out$steer,file.path(dir_res,"17_MR_TF_OA_multi_steiger.csv"))

cat("\n===== PRIMARY MR (lead SNP / TF, IVW = Wald ratio) =====\n")
print(lead_out$main[method=="IVW (RE)", .(tf, outcome, tissue, n_snp,
      est=signif(estimate,3), se=signif(se,3), p=signif(p,3), p_het=signif(p_het,3))])

cat("\n===== SENSITIVITY MR (all cis-eQTLs, IVW) =====\n")
print(multi_out$main[method=="IVW (RE)", .(tf, outcome, tissue, n_snp,
      est=signif(estimate,3), se=signif(se,3), p=signif(p,3), p_het=signif(p_het,3))])

## ---------- 6. figures ----------
make_forest <- function(tab, fname, title){
  ivw <- tab[method=="IVW (RE)" & !is.na(estimate)]
  p <- ggplot(ivw, aes(x=estimate, y=interaction(tf, outcome, sep=" | "))) +
    geom_vline(xintercept=0, linetype="dashed", color="grey40") +
    geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se),
                  orientation="y", height=0.2, color="#2166AC") +
    geom_point(size=2.5, color="#B2182B") +
    facet_wrap(~outcome, scales="free_y") +
    labs(x="MR causal estimate (log OR per 1 SD TF expression), IVW(RE)/Wald",
         y=NULL, title=title) +
    theme_bw(base_size=10) + theme(strip.text=element_text(face="bold"))
  ggsave(fname, p, width=11, height=7, dpi=300)
  p
}
fig64 <- make_forest(lead_out$main,  file.path(dir_fig,"fig64_MR_TF_lead_forest.png"),
                     "Fig 64. MR (primary): lead cis-eQTL per TF -> OA risk")
fig65 <- make_forest(multi_out$main, file.path(dir_fig,"fig65_MR_TF_multi_forest.png"),
                     "Fig 65. MR (sensitivity): all significant cis-eQTLs per TF -> OA risk")

# summary tables for manuscript
summ_lead <- lead_out$main[method=="IVW (RE)", .(tf, outcome, tissue, n_snp,
   est=signif(estimate,3), CI=paste0("[",signif(estimate-1.96*se,2),", ",signif(estimate+1.96*se,2),"]"),
   p=ifelse(p<0.001,"<0.001",signif(p,3)))]
summ_multi <- multi_out$main[method=="IVW (RE)", .(tf, outcome, tissue, n_snp,
   est=signif(estimate,3), CI=paste0("[",signif(estimate-1.96*se,2),", ",signif(estimate+1.96*se,2),"]"),
   p=ifelse(p<0.001,"<0.001",signif(p,3)))]
fwrite(summ_lead,  file.path(dir_res,"17_MR_TF_OA_lead_summary.csv"))
fwrite(summ_multi, file.path(dir_res,"17_MR_TF_OA_multi_summary.csv"))
cat("\nDone. fig64/fig65 saved; CSVs written.\n")
