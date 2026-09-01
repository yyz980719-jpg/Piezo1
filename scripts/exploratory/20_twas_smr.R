#!/usr/bin/env Rscript
# 20_twas_smr.R  -- Phase 2-D: TWAS-lite via SMR (summary-data MR)
#   Exposure: eQTLGen blood cis-eQTL (Z_x per SNP)  [local]
#   Outcome : FinnGen R12 M13 arthrosis (Z_y per SNP) [local]
#   Test    : SMR T = Z_x*Z_y / sqrt(Z_x^2 + Z_y^2) ~ N(0,1) under H0 (no causal effect)
#   Gene-level: top |T| SNP per candidate gene (most associated).
# Note: full genome-wide TWAS (FUSION/S-PrediXcan) needs predictDB weights (not local);
#       this SMR gene-set test is the feasible local alternative.
suppressMessages({ library(data.table) })
source("scripts/_bootstrap.R")
ROOT <- PIEZO1_ROOT
EQ   <- file.path(ROOT, "data/mr/eqtlgen_cis_signif.txt.gz")
GW   <- data_path("mr", "finngen_R12_M13_ARTHROSIS.gz")
OUTC <- file.path(ROOT, "results/20_twas_smr.csv")
OUTF <- file.path(ROOT, "figs/fig69_twas_smr.png")

genes <- c("TCF7L1","GATAD2A","BCL6","ZNF853","PKNOX2","ATF6","PRDM16","ZNF92","PIEZO1")
cat("Loading eQTLGen cis-eQTL for candidate genes...\n"); flush.console()
eq <- fread(EQ, showProgress = FALSE)
eqg <- eq[GeneSymbol %in% genes, .(SNP, GeneSymbol, Zx = Zscore, eAllele = AssessedAllele)]
cat(sprintf("  candidate eQTL SNP-gene pairs: %d\n", nrow(eqg))); flush.console()

cat("Loading OA GWAS and subsetting to candidate SNPs...\n"); flush.console()
gw <- fread(GW, select = c("rsids","ref","alt","beta","sebeta"), showProgress = FALSE)
setnames(gw, "rsids", "SNP")
gws <- gw[SNP %in% unique(eqg$SNP)]
cat(sprintf("  matched in OA GWAS: %d SNPs\n", nrow(gws))); flush.console()

m <- merge(eqg, gws, by = "SNP", all.x = TRUE)
m <- m[!is.na(beta)]
# harmonize effect allele: eQTLGen Z is for eAllele; FinnGen beta is for ALT allele
m[, Zy := beta / sebeta]
m[, Zxh := fifelse(eAllele == alt, Zx, fifelse(eAllele == ref, -Zx, NA_real_))]
m <- m[!is.na(Zxh)]
# SMR per SNP
m[, Tsmr := (Zxh * Zy) / sqrt(Zxh^2 + Zy^2)]
m[, psmr := 2 * pnorm(-abs(Tsmr))]
# gene-level: top |T| SNP
top <- m[, .SD[which.max(abs(Tsmr))], by = GeneSymbol]
out <- top[, .(gene = GeneSymbol, top_snp = SNP, Zx = round(Zxh, 3), Zy = round(Zy, 3),
               T_SMR = round(Tsmr, 3), p_SMR = psmr, eAllele, alt, ref)]
out[, p_fdr := p.adjust(p_SMR, method = "BH")]
setorder(out, p_SMR)
fwrite(out, OUTC)
cat("Wrote", OUTC, "\n"); print(out)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressMessages(library(ggplot2))
  p <- ggplot(out, aes(x = reorder(gene, -abs(T_SMR)), y = T_SMR, fill = -log10(p_SMR))) +
    geom_col(color = "black") +
    scale_fill_gradient(low = "steelblue", high = "red") +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = "dashed", color = "grey40") +
    labs(title = "Candidate-gene TWAS (SMR): blood eQTL -> OA risk (Phase 2-D)",
         subtitle = "eQTLGen cis-eQTL x FinnGen R12 arthrosis; SMR T~N(0,1)",
         x = "gene", y = "SMR test statistic (Z)") +
    theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(OUTF, p, width = 7, height = 5, dpi = 150)
  cat("Wrote", OUTF, "\n")
}
cat("DONE Phase2-D\n")
