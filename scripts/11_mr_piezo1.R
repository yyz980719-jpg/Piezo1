# =====================================================================
# 11_mr_piezo1.R  孟德尔随机化: PIEZO1 表达 -> 骨关节炎(OA)风险
# 暴露: eQTLGen cis-eQTL (全血, N~31,684) 工具变量 6 个 (p<5e-8, r2<0.001)
# 结局: FinnGen R12 M13_ARTHROSIS / _KNEE / _COX (logistic, GRCh38)
# 方法: Wald ratio, IVW(RE), 加权中位数, MR-Egger, Q 异质性, 留一法,
#       近似 Steiger 方向性检验
# 输出: results/11_*.csv, figs/fig43-45_*.png
# =====================================================================

source("scripts/_bootstrap.R")
setwd(data_path("mr"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(20250701)

dir_res <- PIEZO1_RESULTS
dir_fig <- PIEZO1_FIGURES

## ---------- 1. 工具变量 ----------
eq <- fread("eqtlgen_PIEZO1_cis.tsv")
iv <- c("rs56158123","rs137932643","rs61649074","rs149194020","rs74466939","rs77149258")
ivd <- eq[SNP %in% iv]
ivd[, N := as.numeric(NrSamples)]
ivd[, beta_x := Zscore / sqrt(N)]
ivd[, se_x := 1 / sqrt(N)]
ivd[, Fstat := Zscore^2]
cat("工具变量数:", nrow(ivd), "  最小F =", round(min(ivd$Fstat),1), "\n")

## ---------- 2. 结局数据 ----------
outcomes <- list(
  ART = list(file="finngen_M13_ARTHROSIS_chr16region.tsv",  label="OA (overall, M13_ARTHROSIS)",
             ncase=101454, nctrl=315115),
  KNEE= list(file="finngen_M13_ARTHROSIS_KNEE_chr16region.tsv", label="Knee OA (gonarthrosis)",
             ncase=61356,  nctrl=315115),
  COX = list(file="finngen_M13_ARTHTROSIS_COX_chr16region.tsv",  label="Hip OA (coxarthrosis)",
             ncase=30802,  nctrl=315115)
)

## ---------- 3. MR 核心函数 ----------
ivw_re <- function(r, se) {  # 随机效应 IVW (元分析)
  w <- 1/se^2
  k <- length(r); if (k == 1) return(c(est=r, se=se))
  f <- sum(w*r)/sum(w)
  Q <- sum(w*(r-f)^2)
  tau2 <- max(0, (Q - k + 1)/(sum(w) - sum(w^2)/sum(w)))
  wre <- 1/(se^2 + tau2)
  est <- sum(wre*r)/sum(wre); s <- sqrt(1/sum(wre))
  c(est=est, se=s, Q=Q, df=k-1, tau2=tau2)
}
wmedian <- function(r, se, B=5000) {  # 参数化(epsilon) bootstrap, 同时考虑每个比值的不确定性
  w <- 1/se^2; o <- order(r); cw <- cumsum(w[o])/sum(w)
  est <- r[o][which(cw >= 0.5)[1]]
  rb <- replicate(B, { rj <- rnorm(length(r), r, se)
    oo <- order(rj); cww <- cumsum(w[oo])/sum(w); rj[oo][which(cww>=0.5)[1]] })
  c(est=est, se=sd(rb, na.rm=TRUE))
}
mr_egger <- function(bx, by, se_y) {
  fit <- summary(lm(by ~ bx, weights = 1/se_y^2))
  c(slope=coef(fit)[2,1], se_slope=coef(fit)[2,2],
    intercept=coef(fit)[1,1], se_int=coef(fit)[2,2],
    p_int=coef(fit)[1,4])
}

## ---------- 4. 逐结局 MR ----------
# 加载 eQTLGen 等位基因频率 (AlleleB_all) 用于回文位点抢救与频率一致性质控
af_full <- fread("eqtlgen_AF.txt.gz")
af <- af_full[SNP %in% ivd$SNP, .(SNP, AlleleA, AlleleB, AlleleB_all)]
ivd <- merge(ivd, af, by = "SNP", all.x = TRUE)
ivd[, af_assessed := ifelse(AssessedAllele == AlleleB, AlleleB_all, 1 - AlleleB_all)]

snp_res <- list(); main_res <- list(); loo_res <- list(); steer_res <- list()

for (oc in names(outcomes)) {
  o <- outcomes[[oc]]
  fg <- fread(o$file)
  setnames(fg, make.names(names(fg)))  # #chrom -> X.chrom
  m <- merge(ivd, fg, by.x = "SNP", by.y = "rsids")   # 仅保留匹配上的工具变量
  cat("\n==", o$label, "==> 匹配 SNP:", nrow(m), "/", nrow(ivd), "\n")

  # 等位基因和谐化: eQTLGen AssessedAllele vs FinnGen alt
  match_same <- m$AssessedAllele == m$alt & m$OtherAllele == m$ref
  match_swap <- m$AssessedAllele == m$ref & m$OtherAllele == m$alt
  # 频率一致性: 判断 AssessedAllele 在 FinnGen 中对应 alt 还是 ref
  diff_same <- abs(m$af_assessed - m$af_alt)
  diff_swap <- abs(m$af_assessed - (1 - m$af_alt))
  af_consistent <- pmin(diff_same, diff_swap) < 0.08
  keep <- (match_same | match_swap) & af_consistent
  cat("频率不一致剔除:", sum(!af_consistent), " 等位基因完全不符:", sum(!match_same & !match_swap), "\n")
  m <- m[keep]
  is_same <- match_same[keep]
  m$beta_x_h <- ifelse(is_same, m$beta_x, -m$beta_x)
  m[, wald := beta / beta_x_h]
  m[, se_wald := sqrt(sebeta^2/beta_x_h^2 + (beta*se_x)^2/beta_x_h^4)]
  m[, Fstat := Zscore^2]
  m[, outcome := o$label]
  snp_res[[oc]] <- m[, .(outcome, SNP, chr=X.chrom, pos=pos, eff_allele=alt, other_allele=ref,
                        eaf=af_alt, beta_x=beta_x_h, se_x, beta_y=beta, se_y=sebeta,
                        p_y=pval, wald, se_wald, Fstat)]

  k <- nrow(m)
  if (k >= 1) {
    iv <- ivw_re(m$wald, m$se_wald)
    wm <- if (k >= 3) wmedian(m$wald, m$se_wald) else c(est=NA, se=NA)
    eg <- if (k >= 3) mr_egger(m$beta_x_h, m$beta, m$sebeta) else rep(NA_real_, 5)
    main_res[[oc]] <- data.table(
      outcome = o$label, n_snp = k,
      method = c("IVW (RE)", "Weighted median", "MR-Egger slope", "MR-Egger intercept"),
      estimate = c(iv["est"], wm["est"], eg["slope"], eg["intercept"]),
      se = c(iv["se"], wm["se"], eg["se_slope"], eg["se_int"]),
      p = 2*pnorm(abs(c(iv["est"], wm["est"], eg["slope"], eg["intercept"])/
                       c(iv["se"], wm["se"], eg["se_slope"], eg["se_int"])), lower.tail=FALSE),
      Q = c(iv["Q"], NA, NA, NA), p_het = c(pchisq(iv["Q"], iv["df"], lower.tail=FALSE), NA, NA, NA)
    )
    # 留一法
    for (i in 1:k) {
      idx <- setdiff(1:k, i)
      lo <- if (length(idx) >= 1) ivw_re(m$wald[idx], m$se_wald[idx])["est"] else NA
      loo_res[[length(loo_res)+1]] <- data.table(outcome=o$label, excluded=m$SNP[i], est=unname(lo))
    }
    # 近似 Steiger: r2 = Z^2/(Z^2+N)
    r2x <- m$Zscore^2/(m$Zscore^2 + m$N)
    Ny <- o$ncase + o$nctrl; Zy <- m$beta/m$sebeta
    r2y <- Zy^2/(Zy^2 + Ny)
    steer_res[[oc]] <- data.table(outcome=o$label, SNP=m$SNP, r2_exposure=r2x, r2_outcome=r2y,
                                  correct_direction = r2x > r2y)
  }
}

snp_tab  <- rbindlist(snp_res, fill=TRUE)
main_tab <- rbindlist(main_res, fill=TRUE)
loo_tab  <- rbindlist(loo_res)
steer_tab<- rbindlist(steer_res)

fwrite(snp_tab,   file.path(dir_res, "11_MR_snp_level.csv"))
fwrite(main_tab,  file.path(dir_res, "11_MR_main_results.csv"))
fwrite(loo_tab,   file.path(dir_res, "11_MR_leave_one_out.csv"))
fwrite(steer_tab, file.path(dir_res, "11_MR_steiger.csv"))

cat("\n===== 主结果 =====\n")
print(main_tab[, .(outcome, method, n_snp,
                   estimate=signif(estimate,3), se=signif(se,3), p=signif(p,3))])

## ---------- 5. 图 ----------
# fig43: 森林图 (每个 SNP Wald ratio + 各方法汇总)
fs <- rbind(
  snp_tab[, .(outcome, label=paste0("  ", SNP), est=wald, lo=wald-1.96*se_wald, hi=wald+1.96*se_wald, is_method=FALSE)],
  main_tab[, .(outcome, label=method, est=estimate, lo=estimate-1.96*se, hi=estimate+1.96*se, is_method=TRUE)]
)
fs[, outcome := factor(outcome, levels=c("OA (overall, M13_ARTHROSIS)","Knee OA (gonarthrosis)","Hip OA (coxarthrosis)"))]
p <- ggplot(fs, aes(x=est, y=reorder(label, as.numeric(is_method)))) +
  geom_vline(xintercept=0, linetype="dashed", color="grey40") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.25,
                 color=ifelse(fs$is_method, "#B2182B", "grey30")) +
  geom_point(aes(color=is_method, size=is_method)) +
  scale_color_manual(values=c("grey30","#B2182B"), guide="none") +
  scale_size_manual(values=c(1.8,3), guide="none") +
  facet_wrap(~outcome, scales="free_y") +
  labs(x="Wald ratio / causal estimate (log OR per SD PIEZO1 expression)",
       y=NULL, title="Fig 43. MR: PIEZO1 expression -> OA risk") +
  theme_bw(base_size=10) +
  theme(strip.text=element_text(face="bold"))
ggsave(file.path(dir_fig, "fig43_MR_forest.png"), p, width=11, height=8, dpi=300)

# fig44: 留一法
lt <- loo_tab[outcome=="OA (overall, M13_ARTHROSIS)"]
base <- main_tab[outcome=="OA (overall, M13_ARTHROSIS)" & method=="IVW (RE)"]$estimate
p <- ggplot(lt, aes(x=est, y=excluded)) +
  geom_vline(xintercept=base, color="#B2182B", linetype="dashed") +
  geom_point(size=2.5, color="#2166AC") +
  labs(x="IVW estimate (leave-one-out)", y="Excluded SNP",
       title="Fig 44. Leave-one-out sensitivity (overall OA)") +
  theme_bw(base_size=11)
ggsave(file.path(dir_fig, "fig44_MR_leaveoneout.png"), p, width=7, height=4.5, dpi=300)

# fig45: 散点 beta_y vs beta_x
p <- ggplot(snp_tab, aes(x=beta_x, y=beta_y)) +
  geom_errorbar(aes(ymin=beta_y-1.96*se_y, ymax=beta_y+1.96*se_y), width=NA, color="grey60") +
  geom_point(size=2.5, color="#B2182B") +
  geom_hline(yintercept=0, color="grey70") + geom_vline(xintercept=0, color="grey70") +
  facet_wrap(~outcome, scales="free") +
  labs(x="eQTL effect on PIEZO1 expression (SD per allele, eQTLGen)",
       y="Effect on OA risk (log OR, FinnGen R12)",
       title="Fig 45. SNP-level effects: exposure vs outcome") +
  theme_bw(base_size=10) + theme(strip.text=element_text(face="bold"))
ggsave(file.path(dir_fig, "fig45_MR_scatter.png"), p, width=11, height=4, dpi=300)

cat("\nDone. fig43-45 saved.\n")
