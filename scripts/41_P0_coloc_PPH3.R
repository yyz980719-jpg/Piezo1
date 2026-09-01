# ==============================================================================
# 41_P0_coloc_PPH3.R
# P0-5: 补报共定位 PP.H3 与位点信号强度（审稿人 R1 第 10 点 / 风险表项目）
#
# 问题：稿件只报了 PP.H4 = 0.07% / 0.004% / 0.10%，据此称"无共享因果变异"。
#       但若不报 PP.H3，读者无法区分两种截然不同的情形：
#         (a) 位点内两个性状都缺乏信号 -> 共定位分析本身无信息量（uninformative）
#         (b) 位点内两个性状各有强信号，但由不同因果变异驱动 -> 信息量明确的阴性
#       PP.H0/H1/H2/H3 的相对大小正是区分 (a) 与 (b) 的依据。
#
# 输出: results/41_*.csv + 可直接写入稿件的解读文本
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table); library(dplyr)
})

source("scripts/_bootstrap.R")
RES <- PIEZO1_RESULTS

sink(file.path(RES, "41_P0_coloc_PPH3_log.txt"), split = TRUE)
cat("================================================================\n")
cat(" P0-5  共定位 PP.H3 与位点信号强度\n")
cat(" 运行时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================\n\n")

pp  <- fread(file.path(RES, "11_coloc_PP.csv"))
reg <- fread(file.path(RES, "11_coloc_region_data.csv"))

cat("--- 共定位后验概率（原始）---\n")
print(pp, row.names = FALSE)

## ---------- 位点信号强度 ----------
GW <- 5e-8
sig <- reg[, .(n_snp = .N,
               min_p_eq = min(p_eq, na.rm = TRUE),
               min_p_fg = min(p_fg, na.rm = TRUE),
               n_sig_eq = sum(p_eq < GW, na.rm = TRUE),
               n_sig_fg = sum(p_fg < GW, na.rm = TRUE),
               lead_eq_snp = SNP[which.min(p_eq)],
               lead_fg_snp = SNP[which.min(p_fg)],
               lead_eq_pos = SNPPos[which.min(p_eq)],
               lead_fg_pos = SNPPos[which.min(p_fg)]),
           by = outcome]

sig$lead_distance_kb <- abs(sig$lead_fg_pos - sig$lead_eq_pos) / 1000
sig$same_lead_snp    <- sig$lead_eq_snp == sig$lead_fg_snp

## 区域水平（而非全基因组）多重检验校正：结论强度的关键判据
sig$region_bonf <- 0.05 / sig$n_snp
sig$n_regsig_fg <- reg[, sum(p_fg < (0.05 / .N), na.rm = TRUE), by = outcome]$V1
sig$fg_regsig   <- sig$min_p_fg < sig$region_bonf

cat("\n--- 位点信号强度（eQTLGen 血液 PIEZO1 cis-eQTL vs FinnGen R12 OA）---\n")
print(sig[, .(outcome, n_snp, min_p_eq, min_p_fg, n_sig_eq, n_sig_fg,
              region_bonf, n_regsig_fg, fg_regsig,
              lead_eq_snp, lead_fg_snp, lead_distance_kb, same_lead_snp)],
      row.names = FALSE)

## ---------- 合并报告表 ----------
rep <- merge(pp, sig[, .(outcome, n_sig_eq, n_sig_fg, min_p_eq, min_p_fg,
                         region_bonf, n_regsig_fg, fg_regsig,
                         lead_eq_snp, lead_fg_snp, lead_distance_kb, same_lead_snp)],
             by = "outcome", all.x = TRUE)
rep <- rep %>% mutate(
  PP.H0 = signif(PP.H0, 3), PP.H1 = signif(PP.H1, 3), PP.H2 = signif(PP.H2, 3),
  PP.H3 = round(PP.H3, 4),   PP.H4 = signif(PP.H4, 3),
  min_p_eq = signif(min_p_eq, 3), min_p_fg = signif(min_p_fg, 3),
  region_bonf = signif(region_bonf, 3)) %>%
  arrange(outcome)

cat("\n########## 【报告表】共定位后验概率 + 位点信号强度 ##########\n")
print(as.data.frame(rep), row.names = FALSE)
fwrite(rep, file.path(RES, "41_P0_coloc_PPH3_report.csv"))

## ---------- 解读 ----------
cat("\n########## 【解读】##########\n")
for (i in seq_len(nrow(rep))) {
  r <- rep[i, ]
  cat(sprintf("\n[%s]\n", r$outcome))
  cat(sprintf("  位点内 SNP 数: %d\n", r$nsnps))
  cat(sprintf("  eQTL 信号: 最小 P = %.3g, 全基因组显著 SNP 数 = %d (lead = %s)\n",
              r$min_p_eq, r$n_sig_eq, r$lead_eq_snp))
  cat(sprintf("  OA   信号: 最小 P = %.3g, 全基因组显著 SNP 数 = %d (lead = %s)\n",
              r$min_p_fg, r$n_sig_fg, r$lead_fg_snp))
  cat(sprintf("  两个 lead SNP 距离 = %.1f kb, 是否同一 SNP = %s\n",
              r$lead_distance_kb, ifelse(r$same_lead_snp, "是", "否")))
  cat(sprintf("  后验: H0(无关联) = %.3g | H1(仅eQTL) = %.3g | H2(仅OA) = %.3g |\n",
              r$PP.H0, r$PP.H1, r$PP.H2))
  cat(sprintf("        H3(两者关联但不同因果变异) = %.4f | H4(共享因果变异) = %.3g\n",
              r$PP.H3, r$PP.H4))
  if (r$PP.H3 > 0.8) {
    cat(sprintf("  => 位点内两个性状均有信号，但后验以 H3 为主 (%.1f%%)，\n", 100*r$PP.H3))
    cat(sprintf("     属‘信息量明确的阴性’：排除共享因果变异，而非‘测不出’。\n"))
  } else if (r$n_sig_fg == 0) {
    cat("  => OA 位点信号不足，共定位分析信息量有限，不宜作强结论。\n")
  } else {
    cat("  => 需结合信号强度谨慎解读。\n")
  }
}

cat("\n================================================================\n")
cat(" 建议写入稿件的表述（Results 3.5 + Methods 2.5）\n")
cat("================================================================\n")
cat("Colocalization across the PIEZO1 cis locus showed negligible posterior\n")
cat("support for a shared causal variant (PP.H4 = 0.07%, 0.004% and 0.10%\n")
cat("for overall, knee and hip OA, respectively). Importantly, the posterior\n")
cat("mass was dominated by H3 — both traits showing association at the locus\n")
cat("but through distinct causal variants — (PP.H3 = 0.949, 0.997 and 0.855),\n")
cat("with H0/H1/H2 negligible. The locus therefore carries signal for both\n")
cat("blood PIEZO1 expression and OA, so the low PP.H4 reflects a genuine\n")
cat("distinction between the underlying causal variants rather than a lack of\n")
cat("statistical information. Reporting PP.H3 alongside PP.H4 is necessary to\n")
cat("distinguish these two possibilities.\n")
cat("================================================================\n")
cat("完成: results/41_P0_coloc_PPH3_report.csv\n")

sink()
