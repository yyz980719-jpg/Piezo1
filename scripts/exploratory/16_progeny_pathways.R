# =====================================================================
# 16_progeny_pathways.R  —  Upstream signalling pathways (PROGENy)
# Mechanistic upstream of PIEZO1: which signalling pathways drive the
# OA / PIEZO1 program in bulk cartilage (subchondral bone) RNA-seq.
# Input : data/GSE51588_series_matrix.txt.gz (OA 40 / Normal 10)
# Method: progeny() pathway activities per sample; OA vs Normal (limma);
#         Spearman correlation of each pathway with PIEZO1 expression
# Output: figs/fig61_*.png, results/16_*.csv
# =====================================================================
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(data.table)
  library(progeny); library(ggplot2); library(reshape2); library(RColorBrewer)
})

source("scripts/_bootstrap.R")

## ---- 0. load GSE51588 (same as 07/12) ----
g   <- getGEO(filename = "data/GSE51588_series_matrix.txt.gz", getGPL = FALSE)
e   <- as.matrix(exprs(g)); pd <- pData(g)
ttl <- as.character(pd$title)
dis <- ifelse(grepl("^OA-", ttl), "OA", "Normal")
reg <- ifelse(grepl("-MT-", ttl), "Medial", "Lateral")
don <- sub("^[^-]+-[^-]+-", "", ttl)
uid <- paste(dis, don, sep = "_")
if (max(e, na.rm = TRUE) > 100) e <- log2(e + 1)
ann <- fread("data/GPL13497_annot_full.tsv", header = FALSE, skip = 1,
             sep = "\t", colClasses = "character",
             col.names = c("probe", "symbol"))
ann <- ann[!is.na(symbol) & nzchar(symbol) &
           grepl("^[A-Za-z][A-Za-z0-9@\\.\\-]*$", symbol)]
common <- intersect(rownames(e), ann$probe)
eg <- e[common, , drop = FALSE]
gl <- ann$symbol[match(common, ann$probe)]
me <- rowMeans(eg, na.rm = TRUE); o <- order(-me)
eg <- eg[o, , drop = FALSE]; gl <- gl[o]
kd <- !duplicated(gl)
emat <- eg[kd, , drop = FALSE]; rownames(emat) <- gl[kd]
storage.mode(emat) <- "numeric"
pz <- emat["PIEZO1", ]
cat(sprintf("[load] samples=%d (OA %d / Normal %d) genes=%d\n",
            ncol(emat), sum(dis=="OA"), sum(dis=="Normal"), nrow(emat)))

## ---- 1. PROGENy pathway activities ----
pa <- progeny(emat, scale = TRUE, organism = "Human", top = 100)
cat(sprintf("[progeny] activity matrix dims: %d x %d (rows=samples? %s)\n",
            nrow(pa), ncol(pa), ifelse("PIEZO1" %in% rownames(pa), "genes!", "samples")))
# progeny returns pathways (rows) x samples (cols) actually -> ensure samples in rows
if (nrow(pa) > ncol(pa) && !("OA" %in% rownames(pa))) pa <- t(pa)
pa <- as.data.frame(t(pa))            # samples (rows) x pathways (cols)
pa$Sample  <- colnames(emat)
pa$Group   <- dis
pa$Region  <- reg
fwrite(pa, "results/16_pathway_activity.csv")

## ---- 2. OA vs Normal (limma on pathway activities) ----
dsg <- model.matrix(~ factor(dis, levels = c("Normal","OA")))
fit <- eBayes(lmFit(t(as.matrix(pa[, 1:14])), dsg))
tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt$pathway <- rownames(tt)
tt$Group    <- "OA vs Normal"
fwrite(tt[, c("pathway","logFC","AveExpr","t","P.Value","adj.P.Val")],
        "results/16_pathway_OAvsNormal.csv")
cat(sprintf("[OA vs Normal] pathways tested=%d | FDR<0.05: %d\n",
            nrow(tt), sum(tt$adj.P.Val < 0.05)))

## ---- 3. correlation of each pathway with PIEZO1 ----
rho <- sapply(1:14, function(j) suppressWarnings(cor(pa[, j], pz, method = "spearman")))
pcorr <- data.frame(pathway = colnames(pa)[1:14], rho = as.numeric(rho),
                    stringsAsFactors = FALSE)
pcorr$p <- sapply(1:14, function(j)
  suppressWarnings(cor.test(pa[, j], pz, method = "spearman")$p.value))
pcorr$fdr <- p.adjust(pcorr$p, "BH")
pcorr <- pcorr[order(-abs(pcorr$rho)), ]
fwrite(pcorr, "results/16_pathway_PIEZO1_corr.csv")
cat(sprintf("[PIEZO1 corr] pathways |rho|>=0.3 & FDR<0.05: %d\n",
            sum(abs(pcorr$rho)>=0.3 & pcorr$fdr<0.05)))

## ---- FIG 61: OA vs Normal pathway activity (logFC) ----
tt$col <- ifelse(tt$adj.P.Val < 0.05,
                 ifelse(tt$logFC > 0, "#B2182B", "#2166AC"), "grey70")
tt$pathway <- factor(tt$pathway, levels = tt$pathway[order(tt$logFC)])
p61 <- ggplot(tt, aes(x = pathway, y = logFC, fill = col)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  coord_flip() + theme_minimal(base_size = 10) +
  labs(title = "Fig S13A. Upstream signalling pathways altered in OA (PROGENy)",
       subtitle = "Bulk cartilage (subchondral bone); OA vs Normal, limma; red=up, blue=down",
       x = NULL, y = "log2 fold change (OA vs Normal)") +
  theme(legend.position = "none")
ggsave("figs/fig61_progeny_OAvsNormal.png", p61, width = 7, height = 6, dpi = 300)

## ---- FIG 62: pathway activity vs PIEZO1 (Spearman) ----
pcorr$pathway <- factor(pcorr$pathway, levels = pcorr$pathway[order(pcorr$rho)])
p62 <- ggplot(pcorr, aes(x = pathway, y = rho, fill = fdr < 0.05)) +
  geom_col(width = 0.8, color = "white") +
  scale_fill_manual(values = c("TRUE"="#B2182B","FALSE"="grey60"),
                    name = "FDR<0.05") +
  coord_flip() + theme_minimal(base_size = 10) +
  labs(title = "Fig S13B. Pathway activity correlated with PIEZO1 expression",
       subtitle = "Spearman rho across OA+Normal cartilage samples",
       x = NULL, y = "Spearman rho (pathway activity vs PIEZO1)") +
  theme(legend.position = "bottom")
ggsave("figs/fig62_progeny_PIEZO1_corr.png", p62, width = 7, height = 6, dpi = 300)

## ---- FIG 63: heatmap of mean pathway activity by group ----
pa.long <- melt(pa, id.vars = c("Sample","Group","Region"),
                measure.vars = colnames(pa)[1:14],
                variable.name = "pathway", value.name = "act")
setDT(pa.long)
agg <- pa.long[, .(mean_act = mean(act)), by = .(pathway, Group, Region)]
agg$grp <- paste(agg$Group, agg$Region)
wid <- dcast(agg, pathway ~ grp, value.var = "mean_act")
rownames(wid) <- wid$pathway; wid$pathway <- NULL
wid <- as.matrix(wid)
wid <- wid[order(apply(wid, 1, mean)), ]
wid_s <- t(scale(t(wid)))   # z-score per pathway
pdf(NULL)  # ensure device
p63 <- ggplot(melt(wid_s), aes(x = Var2, y = Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "z-score") +
  theme_minimal(base_size = 9) +
  labs(title = "Fig S13C. Mean pathway activity by OA status and region",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figs/fig63_progeny_group_heatmap.png", p63, width = 7, height = 6, dpi = 300)

cat("\n[16] Done. figs fig61-63 + results/16_*.csv written.\n")
