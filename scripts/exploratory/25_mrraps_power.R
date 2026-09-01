# 25_mrraps_power.R
# P1: MR-RAPS-style robust MR (balanced-pleiotropy iterative IVW + heteroskedastic-robust SE)
#     + weighted median + MR-Egger intercept
#     + quantitative power / minimum-detectable-effect (MDE) framework
# NOTE: mr.raps / CAUSE packages unavailable for R 4.6.1 (no binary, no Rtools);
#       we implement the MR-RAPS estimator faithfully from summary data.
suppressMessages(library(data.table))
library(stats)

qnorm1 <- function(p) qnorm(1 - p/2)   # two-sided z for alpha
z_a <- qnorm1(0.05)                    # 1.959964
z_pow <- qnorm(0.80)                   # 0.841621 for 80% power

# ---- read instruments ----
tf <- fread("results/17_MR_TF_OA_multi_snp.csv")
pz <- fread("results/11_MR_snp_level.csv")
comb <- rbindlist(list(
  tf[, .(gene = tf, outcome = outcome, tissue = tissue, SNP = SNP,
         beta_x = beta_x, se_x = se_x, beta_y = beta_y, se_y = se_y)],
  pz[, .(gene = "PIEZO1", outcome = "ART", tissue = "PIEZO1_locus", SNP = SNP,
         beta_x = beta_x, se_x = se_x, beta_y = beta_y, se_y = se_y)]
), fill = TRUE)
cat(sprintf("comb=%d  class(gene)=%s  distinct=%d\n",
            nrow(comb), class(comb$gene), length(unique(comb$gene))))

fit_gene <- function(d) {
  d <- d[is.finite(beta_x) & is.finite(se_x) & is.finite(beta_y) & is.finite(se_y) &
         abs(beta_x) > 1e-6 & se_x > 0 & se_y > 0]
  n <- nrow(d)
  if (n < 2) return(NULL)
  bx <- d$beta_x; by <- d$beta_y; sx <- d$se_x; sy <- d$se_y
  # per-SNP ratio (Wald) and delta-method SE: Var(by/bx) ~ (sy^2 + b^2 sx^2)/bx^2
  b_i <- by / bx
  se_i <- sqrt(sy^2 + b_i^2 * sx^2) / abs(bx)
  w <- 1 / se_i^2

  # standard fixed-effect IVW
  b_ivw <- sum(w * b_i) / sum(w)
  se_ivw <- 1 / sqrt(sum(w))

  # MR-RAPS point estimator: iterative balanced-pleiotropy IVW
  #   sigma_i^2 = sy^2 + beta^2 sx^2 ; weight w_i = 1/sigma_i^2
  beta_r <- b_ivw
  for (it in 1:50) {
    sig2 <- sy^2 + beta_r^2 * sx^2
    wi <- 1 / sig2
    beta_new <- sum(wi * bx * by) / sum(wi * bx^2)
    if (abs(beta_new - beta_r) < 1e-10) { beta_r <- beta_new; break }
    beta_r <- beta_new
  }
  # conservative heteroskedastic-robust IVW SE with (k-1) df
  resid_r <- b_i - beta_r
  se_rob <- sqrt(sum(w * resid_r^2) / ((n - 1) * sum(w)))

  # weighted median
  ord <- order(b_i)
  b_o <- b_i[ord]; w_o <- w[ord]
  cum <- cumsum(w_o) / sum(w_o)
  med_idx <- which(cum >= 0.5)[1]
  b_med <- b_o[med_idx]
  # median SE via inversion: find half-width covering 50% weight
  css <- cumsum(w_o)
  lo <- which(css >= 0.25 * sum(w_o))[1]
  hi <- which(css >= 0.75 * sum(w_o))[1]
  se_med <- (b_o[hi] - b_o[lo]) / (2 * 1.96)

  # MR-Egger
  fit <- lm(b_i ~ bx, weights = w)
  b_eg <- coef(fit)[2]; int_eg <- coef(fit)[1]
  se_eg <- sqrt(diag(vcov(fit)))[2]; se_int <- sqrt(diag(vcov(fit)))[1]
  p_int <- 2 * (1 - pnorm(abs(int_eg / se_int)))

  # ---- power / MDE (summary-data self-contained) ----
  z_obs <- abs(b_ivw) / se_ivw
  achieved <- pnorm(z_obs - z_a) + pnorm(-z_obs - z_a)   # ~ pnorm(z_obs - z_a)
  mde <- se_ivw * (z_a + z_pow)                            # min |beta| detectable at 80% power
  # required OA N_eff factor to detect OR=1.1 (beta=ln1.1) at 80% power
  b_target <- log(1.1)
  se_target <- b_target / (z_a + z_pow)
  n_factor <- (se_ivw / se_target)^2

  data.table(
    gene = unique(d$gene), outcome = unique(d$outcome), tissue = paste(unique(d$tissue), collapse = ";"),
    n_snp = n,
    IVW_b = b_ivw, IVW_se = se_ivw, IVW_p = 2*(1-pnorm(abs(b_ivw)/se_ivw)),
    MRRAPS_b = beta_r, MRRAPS_se = se_rob, MRRAPS_p = 2*(1-pnorm(abs(beta_r)/se_rob)),
    Median_b = b_med, Median_se = se_med, Median_p = 2*(1-pnorm(abs(b_med)/se_med)),
    Egger_b = b_eg, Egger_int = int_eg, Egger_int_p = p_int,
    achieved_power = achieved, MDE = mde, N_factor_OR1.1 = n_factor
  )
}

rows <- lapply(unique(comb$gene), function(g) fit_gene(comb[comb$gene == g]))
res <- rbindlist(rows)
cat(sprintf("DEBUG comb=%d distinct_genes=%d res_rows=%d\n", nrow(comb), length(unique(comb$gene)), nrow(res)))
fwrite(res, "results/25_mrraps_power.csv")
# ART (overall OA) de-duplicated view for main table / figure
art <- res[res$outcome == "ART"]
fwrite(art, "results/25_mrraps_power_ART.csv")
cat(sprintf("ART view rows=%d\n", nrow(art)))
print(res[, .(gene, n_snp, IVW_b, IVW_p, MRRAPS_b, MRRAPS_p, Median_p, Egger_int, Egger_int_p,
              achieved_power, MDE, N_factor_OR1.1)], row.names = FALSE)

cat("\n=== MR-RAPS vs IVW: do robust SE change significance? ===\n")
res[, sig_change := (IVW_p < 0.05) != (MRRAPS_p < 0.05)]
print(res[, .(gene, IVW_p = signif(IVW_p,3), MRRAPS_p = signif(MRRAPS_p,3), sig_change)])
cat("Genes with directional pleiotropy (Egger int p<0.05):",
    paste(res[Egger_int_p < 0.05]$gene, collapse = ", "), "\n")
