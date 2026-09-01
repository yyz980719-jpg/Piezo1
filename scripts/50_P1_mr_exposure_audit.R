# =====================================================================
# 50_P1_mr_exposure_audit.R
# 模拟审稿后优化 · T1.1 (MR exposure effect-scale audit) + T1.2 (MR diagnostics)
#
# 目标
#   1) 审计 eQTLGen exposure 效应量的来源与量纲 (beta = Z/sqrt(N), SE = 1/sqrt(N))
#   2) 暴露原稿 Fstat = Z^2 的循环论证问题，改用 MAF-based F
#   3) 补齐 MR diagnostics: 逐 SNP F/R2、Cochran Q、I2、MR-Egger intercept、
#      Steiger-type 方向性检验、IVW 权重构成、FinnGen case/control 数
#   4) 证明 IVW 推断对 exposure 量纲缩放不变，并对 power/MDE 做量纲误差敏感性
#
# 输出: results/50_*.csv
# =====================================================================

source("scripts/_bootstrap.R")
suppressPackageStartupMessages({ library(data.table) })
set.seed(20250701)

dir_res <- PIEZO1_RESULTS
dir_mr  <- data_path("mr")
sink(result_path("50_P1_mr_exposure_audit_log.txt"), split = TRUE)
cat("=== 50_P1 MR exposure effect-scale audit ===\n")
cat("R:", R.version.string, "| data.table:", as.character(packageVersion("data.table")), "\n\n")

## ---------------------------------------------------------------------
## 0. 结局表型定义 (FinnGen R12 manifest)
## ---------------------------------------------------------------------
out_meta <- data.table(
  outcome = c("OA (overall, M13_ARTHROSIS)", "Knee OA (gonarthrosis)", "Hip OA (coxarthrosis)"),
  phenocode = c("M13_ARTHROSIS", "M13_ARTHROSIS_KNEE", "M13_ARTHTROSIS_COX"),
  n_case = c(101454, 61356, 30802),
  n_ctrl = c(315115, 315115, 315115)
)
out_meta[, n_total := n_case + n_ctrl]
out_meta[, prev := n_case / n_total]

## ---------------------------------------------------------------------
## 1. Exposure: eQTLGen PIEZO1 cis-eQTL —— 字段可用性审计
## ---------------------------------------------------------------------
eq <- fread(file.path(dir_mr, "eqtlgen_PIEZO1_cis.tsv"))
cat("eQTLGen PIEZO1 cis 记录数:", nrow(eq), "\n")
cat("可用字段:", paste(names(eq), collapse = ", "), "\n")
cat("是否包含 beta/SE/EAF:",
    ifelse(all(c("beta","se","EAF") %in% names(eq)), "是", "否 (仅 Z score + N)"), "\n\n")

iv <- c("rs56158123","rs137932643","rs61649074","rs149194020","rs74466939","rs77149258")
ivd <- eq[SNP %in% iv, .(SNP, AssessedAllele, OtherAllele, Zscore, NrSamples, Pvalue)]
ivd[, N := as.numeric(NrSamples)]
ivd[, beta_x := Zscore / sqrt(N)]
ivd[, se_x  := 1 / sqrt(N)]

## 原稿做法: F = Z^2 (循环论证 —— 由重构本身导出, 不提供独立信息)
ivd[, F_naive := Zscore^2]

## ---------------------------------------------------------------------
## 2. Allele frequency -> 非循环论证的 R2 与 F
##    R2 = 2*MAF*(1-MAF)*beta^2   (表型已逆正态标准化, SD = 1)
##    F  = R2 * (N - 2) / (1 - R2)
## ---------------------------------------------------------------------
af <- fread(file.path(dir_mr, "eqtlgen_AF_instruments.tsv"))
af[, freq_A := 1 - AlleleB_all]
af_long <- rbind(
  af[, .(SNP, allele = AlleleA, freq = freq_A)],
  af[, .(SNP, allele = AlleleB, freq = AlleleB_all)]
)
ivd <- merge(ivd, af_long, by.x = c("SNP","AssessedAllele"), by.y = c("SNP","allele"), all.x = TRUE)
setnames(ivd, "freq", "maf_ea")
ivd[, R2_x := 2 * maf_ea * (1 - maf_ea) * beta_x^2]
ivd[, F_maf := R2_x * (N - 2) / (1 - R2_x)]

cat("=== 工具变量与效应量重构 ===\n")
print(ivd[order(-abs(Zscore)), .(SNP, AssessedAllele, OtherAllele, Zscore, N,
                                 beta_x = round(beta_x, 5), se_x = round(se_x, 5),
                                 maf_ea = round(maf_ea, 4),
                                 F_naive_Z2 = round(F_naive, 1),
                                 R2_x = round(R2_x, 6), F_maf = round(F_maf, 2))])
cat("\n>> 关键审计结论: F_naive = Z^2 恒等于重构本身 (F = (Z/sqrtN)^2 / (1/N)), 不提供独立的\n",
    "   工具变量强度信息。改用 MAF-based F 后, 仅 rs56158123 满足 F > 10。\n\n")

## ---------------------------------------------------------------------
## 3. 结局与 SNP 层面数据
## ---------------------------------------------------------------------
snp <- fread(file.path(dir_res, "11_MR_snp_level.csv"))
snp <- merge(snp, ivd[, .(SNP, maf_ea, beta_x2 = beta_x, se_x2 = se_x, F_maf, R2_x, N)],
             by = "SNP", all.x = TRUE)
stopifnot(!any(is.na(snp$maf_ea)))

## ---------------------------------------------------------------------
## 4. MR 重算: IVW(FE/RE)、加权中位数、MR-Egger、Q、I2、留一法、权重构成
## ---------------------------------------------------------------------
ivw_fe <- function(r, se) c(est = sum(r/se^2)/sum(1/se^2), se = sqrt(1/sum(1/se^2)))
ivw_re <- function(r, se) {
  w <- 1/se^2; k <- length(r)
  f <- sum(w*r)/sum(w); Q <- sum(w*(r-f)^2)
  tau2 <- max(0, (Q - k + 1)/(sum(w) - sum(w^2)/sum(w)))
  wre <- 1/(se^2 + tau2)
  c(est = sum(wre*r)/sum(wre), se = sqrt(1/sum(wre)), Q = Q, df = k-1, tau2 = tau2)
}
wmedian <- function(r, se, B = 5000) {
  w <- 1/se^2; o <- order(r); cw <- cumsum(w[o])/sum(w)
  est <- r[o][which(cw >= 0.5)[1]]
  rb <- replicate(B, { rj <- rnorm(length(r), r, se)
    oo <- order(rj); cww <- cumsum(w[oo])/sum(w); rj[oo][which(cww >= 0.5)[1]] })
  c(est = est, se = sd(rb, na.rm = TRUE))
}

mr_rows <- list(); het_rows <- list(); wrows <- list(); loo_rows <- list()

for (i in seq_len(nrow(out_meta))) {
  oc <- out_meta$outcome[i]; d <- snp[outcome == oc]
  d[, wald := beta_y / beta_x]
  d[, se_wald := abs(se_y / beta_x)]
  d[, w_ivw := beta_x^2 / se_y^2]

  fe <- ivw_fe(d$wald, d$se_wald); re <- ivw_re(d$wald, d$se_wald)
  wm <- wmedian(d$wald, d$se_wald)
  eg <- summary(lm(d$beta_y ~ d$beta_x, weights = 1/d$se_y^2))
  Q <- re[["Q"]]; dfq <- re[["df"]]
  I2 <- max(0, (Q - dfq) / Q) * 100

  mr_rows[[length(mr_rows) + 1]] <- data.table(
    outcome = oc, n_snp = nrow(d),
    ivw_fe_beta = fe[["est"]], ivw_fe_se = fe[["se"]],
    ivw_fe_p = 2*pnorm(-abs(fe[["est"]]/fe[["se"]])),
    ivw_re_beta = re[["est"]], ivw_re_se = re[["se"]],
    ivw_re_p = 2*pnorm(-abs(re[["est"]]/re[["se"]])),
    wm_beta = wm[["est"]], wm_se = wm[["se"]], wm_p = 2*pnorm(-abs(wm[["est"]]/wm[["se"]])),
    egger_slope = coef(eg)[2,1], egger_slope_se = coef(eg)[2,2],
    egger_slope_p = coef(eg)[2,4],
    egger_intercept = coef(eg)[1,1], egger_intercept_se = coef(eg)[1,2],
    egger_intercept_p = coef(eg)[1,4])
  het_rows[[length(het_rows) + 1]] <- data.table(
    outcome = oc, Q = Q, df = dfq, p_het = pchisq(Q, dfq, lower.tail = FALSE),
    I2_pct = I2, tau2 = re[["tau2"]])

  wrows[[length(wrows) + 1]] <- data.table(
    outcome = oc, SNP = d$SNP, beta_x = d$beta_x, se_y = d$se_y,
    w_ivw = d$w_ivw, weight_pct = 100*d$w_ivw/sum(d$w_ivw), F_maf = d$F_maf)

  for (j in seq_len(nrow(d))) {
    dd <- d[-j]
    r <- ivw_re(dd$wald, dd$se_wald)
    loo_rows[[length(loo_rows) + 1]] <- data.table(
      outcome = oc, excluded = d$SNP[j], ivw_re_beta = r[["est"]],
      ivw_re_se = r[["se"]], ivw_re_p = 2*pnorm(-abs(r[["est"]]/r[["se"]])))
  }
}
mr  <- rbindlist(mr_rows);  het <- rbindlist(het_rows)
wts <- rbindlist(wrows);    loo <- rbindlist(loo_rows)

cat("=== MR 主结果 (重算) ===\n")
print(mr[, .(outcome, n_snp, ivw_re_beta = round(ivw_re_beta, 4), ivw_re_se = round(ivw_re_se, 4),
             ivw_re_p = signif(ivw_re_p, 3), wm_beta = round(wm_beta, 4), wm_p = signif(wm_p, 3),
             egger_slope = round(egger_slope, 4), egger_intercept_p = signif(egger_intercept_p, 3))])
cat("\n=== 异质性 ===\n")
print(het[, .(outcome, Q = round(Q, 3), df, p_het = round(p_het, 4), I2_pct = round(I2_pct, 1))])
cat("\n=== IVW 权重构成 (识别 MR 是否实质为单 SNP) ===\n")
print(wts[order(outcome, -weight_pct), .(outcome, SNP, F_maf = round(F_maf, 2), weight_pct = round(weight_pct, 2))])
cat("\n>> 关键审计结论: IVW 权重高度集中于 rs56158123 (唯一 F>10 的工具变量)。\n")
cat("   MR 结果实质上近似单 SNP Wald ratio, 必须在稿件中说明。\n\n")

## ---------------------------------------------------------------------
## 5. Steiger-type 方向性检验 (两样本相关系数 Fisher-z 比较)
##    r_x = beta_x * sqrt(2*MAF*(1-MAF))            (标准化连续暴露)
##    r_y = beta_logit*P(1-P) * sqrt(2*MAF*(1-MAF)) / sqrt(P(1-P))   (二分类结局)
## ---------------------------------------------------------------------
st_rows <- list()
for (i in seq_len(nrow(out_meta))) {
  oc <- out_meta$outcome[i]; Pv <- out_meta$prev[i]; Nout <- out_meta$n_total[i]
  d <- snp[outcome == oc]
  sd_g <- sqrt(2*d$maf_ea*(1-d$maf_ea))
  r_x <- d$beta_x * sd_g
  r_y <- (d$beta_y * Pv*(1-Pv)) * sd_g / sqrt(Pv*(1-Pv))
  zf  <- (atanh(r_x) - atanh(r_y)) / sqrt(1/(d$N - 3) + 1/(Nout - 3))
  st_rows[[length(st_rows)+1]] <- data.table(
    outcome = oc, SNP = d$SNP,
    r2_exposure = r_x^2, r2_outcome = r_y^2,
    steiger_z = zf, steiger_p = 2*pnorm(-abs(zf)),
    direction_correct = r_x^2 > r_y^2)
}
steiger <- rbindlist(st_rows)
cat("=== Steiger-type 方向性检验 ===\n")
print(steiger[, .(outcome, SNP, r2_exposure = signif(r2_exposure, 3), r2_outcome = signif(r2_outcome, 3),
                  steiger_z = round(steiger_z, 2), steiger_p = signif(steiger_p, 3), direction_correct)])
cat("\n")

## ---------------------------------------------------------------------
## 6. 量纲不变性证明 + power/MDE 量纲敏感性
##    若真实 beta = c * beta_reconstructed, 则 Wald ratio 缩小 c 倍, SE 同比例缩小,
##    IVW z 统计量不变 -> 全部统计推断 (P 值、Q、Egger intercept P) 不受影响。
##    仅 SD 单位下的效应量与 MDE 受 c 影响。
## ---------------------------------------------------------------------
cat("=== 量纲不变性数值验证 (c = exposure 缩放因子) ===\n")
inv <- rbindlist(lapply(c(0.6, 0.8, 1.0, 1.25, 1.67), function(cf) {
  rbindlist(lapply(seq_len(nrow(out_meta)), function(i) {
    oc <- out_meta$outcome[i]; d <- snp[outcome == oc]
    w  <- (d$beta_y) / (cf*d$beta_x); s <- abs(d$se_y/(cf*d$beta_x))
    fe <- ivw_fe(w, s)
    data.table(c_scale = cf, outcome = oc,
               ivw_beta = fe[["est"]], ivw_se = fe[["se"]],
               z = fe[["est"]]/fe[["se"]], p = 2*pnorm(-abs(fe[["est"]]/fe[["se"]])))
  }))
}))
print(inv[, .(c_scale, outcome, ivw_beta = round(ivw_beta, 4), ivw_se = round(ivw_se, 4),
              z = round(z, 4), p = round(p, 4))])
cat("\n>> z 与 P 在任意缩放因子下完全相同: MR 的‘无证据’结论对 exposure 量纲假设不敏感。\n\n")

## power / MDE (对 OR 的检验, alpha = 0.05, 双侧)
mde <- function(se, power = 0.80) (qnorm(1-0.05/2) + qnorm(power)) * se
pow <- function(se, beta) pnorm(-qnorm(1-0.05/2) + abs(beta)/se) + pnorm(-qnorm(1-0.05/2) - abs(beta)/se)
pw <- rbindlist(lapply(seq_len(nrow(out_meta)), function(i) {
  oc <- out_meta$outcome[i]; d <- snp[outcome == oc]
  w <- d$beta_y/d$beta_x; s <- abs(d$se_y/d$beta_x)
  re <- ivw_re(w, s); fe <- ivw_fe(w, s)
  data.table(outcome = oc, n_case = out_meta$n_case[i], n_ctrl = out_meta$n_ctrl[i],
             ivw_se_re = re[["se"]], ivw_se_fe = fe[["se"]],
             mde_beta_re = mde(re[["se"]]), mde_beta_fe = mde(fe[["se"]]))
}))
pw[, mde_OR_re := exp(mde_beta_re)]
pw[, mde_OR_fe := exp(mde_beta_fe)]
## 主分析用随机效应 SE (存在异质性时更保守, 与原稿口径一致); 固定效应作为敏感性
for (orv in c(1.10, 1.20, 1.50, 2.00)) {
  pw[, paste0("power_OR_", orv) := pow(ivw_se_re, log(orv))]
}
## 量纲误差敏感性: 若真实量纲为重构值的 c 倍, MDE(SD 单位) 变为 mde_beta/c
pw_sens <- rbindlist(lapply(c(0.6, 0.8, 1.0, 1.25, 1.67), function(cf)
  pw[, .(c_scale = cf, outcome, mde_OR = exp(mde_beta_re/cf))]))
cat("=== Power / MDE (主分析: 随机效应 SE) ===\n")
print(pw[, .(outcome, n_case, ivw_se_re = round(ivw_se_re, 4), mde_beta_re = round(mde_beta_re, 3),
             mde_OR_re = round(mde_OR_re, 3), power_OR_1.10 = round(power_OR_1.1, 3),
             power_OR_1.20 = round(power_OR_1.2, 3), power_OR_1.50 = round(power_OR_1.5, 3))])
cat("\n=== 敏感性: 固定效应 SE 下的 MDE ===\n")
print(pw[, .(outcome, ivw_se_fe = round(ivw_se_fe, 4), mde_OR_fe = round(mde_OR_fe, 3))])
cat("\n=== MDE 对 exposure 量纲误差的敏感性 ===\n")
print(dcast(pw_sens, outcome ~ c_scale, value.var = "mde_OR"))
cat("\n")

## ---------------------------------------------------------------------
## 7. 写出审计表
## ---------------------------------------------------------------------
fwrite(ivd, file.path(dir_res, "50_instrument_scale_audit.csv"))
fwrite(mr,  file.path(dir_res, "50_mr_rerun.csv"))
fwrite(het, file.path(dir_res, "50_mr_heterogeneity.csv"))
fwrite(wts, file.path(dir_res, "50_mr_ivw_weights.csv"))
fwrite(loo, file.path(dir_res, "50_mr_leave_one_out.csv"))
fwrite(steiger, file.path(dir_res, "50_mr_steiger.csv"))
fwrite(pw,   file.path(dir_res, "50_mr_power.csv"))
fwrite(pw_sens, file.path(dir_res, "50_mr_power_scale_sensitivity.csv"))
fwrite(inv,  file.path(dir_res, "50_mr_scale_invariance.csv"))
cat("审计表已写出: results/50_*.csv\n")
sink()
