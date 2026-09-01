# =====================================================================
# 11b_coloc_piezo1.R  共定位分析: PIEZO1 cis-eQTL 与 OA GWAS 信号
# eQTLGen 全量 cis 统计 (10,278 SNP, hg19) x FinnGen R12 (GRCh38, 按 rsID 匹配)
# coloc.abf 假设检验: H4 = 共享同一因果变异
# 输出: results/11_coloc_*.csv, figs/fig46-47
# =====================================================================

source("scripts/_bootstrap.R")
setwd(data_path("mr"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(coloc) })

dir_res <- PIEZO1_RESULTS
dir_fig <- PIEZO1_FIGURES

## ---------- 1. eQTLGen 全量 cis ----------
eq <- fread("eqtlgen_PIEZO1_cis_full.tsv")
eq[, N := as.numeric(NrSamples)]
eq[, beta_eq := Zscore / sqrt(N)]
eq[, varbeta_eq := 1 / N]
cat("eQTLGen PIEZO1 cis SNPs:", nrow(eq), "\n")

## ---------- 2. FinnGen ----------
outcomes <- list(
  ART = list(file="finngen_M13_ARTHROSIS_chr16region.tsv",  label="OA (overall)",
             ncase=101454, nctrl=315115),
  KNEE= list(file="finngen_M13_ARTHROSIS_KNEE_chr16region.tsv", label="Knee OA",
             ncase=61356,  nctrl=315115),
  COX = list(file="finngen_M13_ARTHTROSIS_COX_chr16region.tsv",  label="Hip OA",
             ncase=30802,  nctrl=315115)
)

coloc_rows <- list(); regdat <- list()
for (oc in names(outcomes)) {
  o <- outcomes[[oc]]
  fg <- fread(o$file)
  setnames(fg, make.names(names(fg)))
  m <- merge(eq, fg, by.x = "SNP", by.y = "rsids")
  cat("\n==", o$label, "==> 匹配 SNP:", nrow(m), "\n")

  # 等位基因和谐化 (与 MR 一致): 统一到 alt 等位基因方向
  match_same <- m$AssessedAllele == m$alt & m$OtherAllele == m$ref
  match_swap <- m$AssessedAllele == m$ref & m$OtherAllele == m$alt
  pal <- (m$AssessedAllele == "A" & m$OtherAllele == "T") | (m$AssessedAllele == "T" & m$OtherAllele == "A") |
         (m$AssessedAllele == "C" & m$OtherAllele == "G") | (m$AssessedAllele == "G" & m$OtherAllele == "C")
  keep <- (match_same | match_swap)
  cat("回文保留(方向按 alt 统一,敏感性中再剔除):", sum(pal & keep), " 不符剔除:", sum(!keep), "\n")
  m <- m[keep]
  # eQTL beta 统一到 FinnGen alt 方向
  m[, beta_h := ifelse(match_same[keep], m$beta_eq, -m$beta_eq)]

  Ny <- o$ncase + o$nctrl; sfrac <- o$ncase / Ny
  d1 <- list(snp = m$SNP, beta = m$beta_h, varbeta = m$varbeta_eq, N = m$N, sdY = 1, type = "quant")
  d2 <- list(snp = m$SNP, beta = m$beta, varbeta = m$sebeta^2, N = Ny, s = sfrac, type = "cc")
  res <- suppressWarnings(coloc.abf(dataset1 = d1, dataset2 = d2))

  pp <- as.data.table(t(res$summary))
  cat(sprintf("%s: PP.H4 = %.3f | nsnps = %d\n", o$label, pp$PP.H4.abf, pp$nsnps))
  coloc_rows[[oc]] <- data.table(outcome = o$label, nsnps = pp$nsnps,
    PP.H0 = pp$PP.H0.abf, PP.H1 = pp$PP.H1.abf, PP.H2 = pp$PP.H2.abf,
    PP.H3 = pp$PP.H3.abf, PP.H4 = pp$PP.H4.abf)
  m[, outcome := o$label]
  regdat[[oc]] <- m[, .(outcome, SNP, SNPPos, p_eq = Pvalue, p_fg = pval, beta_eq = beta_h, beta_fg = beta)]
}

coloc_tab <- rbindlist(coloc_rows)
fwrite(coloc_tab, file.path(dir_res, "11_coloc_PP.csv"))
regall <- rbindlist(regdat)
fwrite(regall, file.path(dir_res, "11_coloc_region_data.csv"))

## ---------- 3. 图 ----------
# fig46: 区域并列曼哈顿 (总体 OA)
rg <- regall[outcome == "OA (overall)"]
rg[, neglog_p_eq := -log10(p_eq)][, neglog_p_fg := -log10(p_fg)]
gmax <- max(rg$SNPPos); gmin <- min(rg$SNPPos)
eq_top <- rg[which.min(p_eq)]; fg_top <- rg[which.min(p_fg)]
lab <- rbind(
  data.table(SNPPos = eq_top$SNPPos, y = eq_top$neglog_p_eq, lab = paste0(eq_top$SNP, "\n(eQTLGen P=", signif(eq_top$p_eq,2), ")"), panel="A. eQTLGen PIEZO1 cis-eQTL"),
  data.table(SNPPos = fg_top$SNPPos, y = fg_top$neglog_p_fg, lab = paste0(fg_top$SNP, "\n(FinnGen P=", signif(fg_top$p_fg,2), ")"), panel="B. FinnGen overall OA"))
p <- ggplot(rg, aes(x = SNPPos/1e6)) +
  geom_point(aes(y = neglog_p_eq), size = 0.5, color = "#2166AC") +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey50") +
  geom_text(data = lab[panel == "A. eQTLGen PIEZO1 cis-eQTL"], aes(y = y, label = lab), size = 2.8, vjust = 1.2, inherit.aes = TRUE) +
  labs(x = "Chromosome 16 position (Mb, hg19)", y = "-log10(P)",
       title = "Fig 46. PIEZO1 cis-eQTL signal (eQTLGen)") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, "fig46_coloc_region_eqtl.png"), p, width = 8, height = 4, dpi = 300)
p2 <- ggplot(rg, aes(x = SNPPos/1e6)) +
  geom_point(aes(y = neglog_p_fg), size = 0.5, color = "#B2182B") +
  geom_text(data = lab[panel == "B. FinnGen overall OA"], aes(y = y, label = lab), size = 2.8, vjust = 1.2, inherit.aes = TRUE) +
  labs(x = "Chromosome 16 position (Mb, hg19)", y = "-log10(P)",
       title = "Fig 46. OA GWAS signal at PIEZO1 locus (FinnGen)") +
  theme_bw(base_size = 11)
ggsave(file.path(dir_fig, "fig46b_coloc_region_oa.png"), p2, width = 8, height = 4, dpi = 300)

# fig47: locuscompare 风格散点 (-log10 P vs -log10 P)
p <- ggplot(rg, aes(x = neglog_p_eq, y = neglog_p_fg)) +
  geom_point(size = 0.8, alpha = 0.6, color = "#4A7BB7") +
  geom_smooth(method = "loess", se = FALSE, color = "#B2182B", linewidth = 0.8) +
  facet_wrap(~outcome, scales = "free") +
  labs(x = "eQTLGen PIEZO1 cis-eQTL  -log10(P)",
       y = "FinnGen OA GWAS  -log10(P)",
       title = "Fig 47. Locuscompare: signal concordance at PIEZO1 locus") +
  theme_bw(base_size = 10) + theme(strip.text = element_text(face = "bold"))
ggsave(file.path(dir_fig, "fig47_coloc_locuscompare.png"), p, width = 11, height = 4, dpi = 300)

cat("\nColoc done.\n"); print(coloc_tab)
