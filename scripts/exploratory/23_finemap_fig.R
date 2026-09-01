#!/usr/bin/env Rscript
# 23_finemap_fig.R -- locus fine-mapping overview (Fig 9)
# Reads per-locus TSVs from 23_finemap_abf.py, plots -log10(p) vs position,
# highlights the 95% credible set (orange), the GWAS lead SNP (red),
# and the eQTL-instrument center (blue dashed).
library(ggplot2)
library(data.table)

loci <- list(
  PIEZO1 =list(file="results/23_finemap_PIEZO1.tsv", center=88746077),
  TCF7L1 =list(file="results/23_finemap_TCF7L1.tsv", center=85264624),
  GATAD2A=list(file="results/23_finemap_GATAD2A.tsv", center=19529715),
  BCL6   =list(file="results/23_finemap_BCL6.tsv",   center=187983973),
  ZNF853 =list(file="results/23_finemap_ZNF853.tsv", center=6499492),
  PKNOX2 =list(file="results/23_finemap_PKNOX2.tsv", center=125090281),
  ATF6   =list(file="results/23_finemap_ATF6.tsv",   center=161942933)
)
dfs <- list()
for (g in names(loci)) {
  d <- fread(loci[[g]]$file)
  d$gene <- g
  d$center <- loci[[g]]$center
  dfs[[g]] <- d
}
D <- rbindlist(dfs)
D <- D[D$pval > 0 & is.finite(D$pval), ]
D$nlp <- -log10(D$pval)
lead <- D[D$pval == min(pval), .SD[1], by = gene]

p <- ggplot(D, aes(x = pos / 1e6, y = nlp)) +
  geom_point(aes(color = factor(in_credible_set)), size = 0.35, alpha = 0.45) +
  geom_vline(aes(xintercept = center / 1e6), linetype = "dashed", color = "steelblue") +
  geom_point(data = lead, aes(x = pos / 1e6, y = nlp), color = "red", size = 1.6, shape = 17) +
  geom_hline(yintercept = 7.3, linetype = "dotted", color = "#cc0000") +
  facet_wrap(~gene, scales = "free_x", ncol = 3) +
  scale_color_manual(values = c("grey70", "darkorange"), guide = "none") +
  labs(x = "Genomic position (Mb)", y = expression(-log[10](p)),
       title = "Locus fine-mapping of PIEZO1 and six consensus transcription-factor loci",
       subtitle = "FinnGen R12 OA (any-site). Orange = 95% ABF credible set; red triangle = GWAS lead SNP; blue dashed = cis-eQTL instrument center; dotted = genome-wide threshold") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "grey90"),
        plot.title = element_text(size = 11, face = "bold"))

ggsave("figs/fig72_finemap_locus.png", p, width = 10, height = 7.5, dpi = 150)
cat("saved figs/fig72_finemap_locus.png\n")
