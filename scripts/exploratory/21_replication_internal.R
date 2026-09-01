#!/usr/bin/env Rscript
# 21_replication_internal.R  -- Phase 2-A (interim): within-FinnGen R12 internal replication
#   Discovery = overall OA (ART); replications = knee (KNEE) and hip (COX) arthrosis,
#   three partially independent OA definitions from the SAME R12 GWAS round.
#   Uses the B-pipeline multi-SNP IVW(RE) results already on disk (no new download).
suppressMessages({ library(data.table); library(ggplot2) })
source("scripts/_bootstrap.R")
ROOT <- PIEZO1_ROOT
IN  <- file.path(ROOT, "results/17_MR_TF_OA_multi_main.csv")
OUTC<- file.path(ROOT, "results/21_replication_internal.csv")
OUTF<- file.path(ROOT, "figs/fig70_replication_internal.png")

d <- fread(IN)[method == "IVW (RE)"]
d <- d[outcome %in% c("ART","KNEE","COX")]
# wide: estimate & p per outcome
w <- dcast(d, tf ~ outcome, value.var = "estimate")
p <- dcast(d, tf ~ outcome, value.var = "p")
colnames(w) <- ifelse(colnames(w)=="tf","tf", paste0(colnames(w),"_est"))
colnames(p) <- ifelse(colnames(p)=="tf","tf", paste0(colnames(p),"_p"))
m <- merge(w, p, by="tf")
setnames(m, c("ART_est","KNEE_est","COX_est","ART_p","KNEE_p","COX_p"),
            c("overall_beta","knee_beta","hip_beta","overall_p","knee_p","hip_p"))

# direction consistency (all non-NA same sign)
m[, dir_consist := sign(overall_beta)==sign(knee_beta) & sign(overall_beta)==sign(hip_beta)]
m[, n_sig05 := rowSums(cbind(overall_p<0.05, knee_p<0.05, hip_p<0.05), na.rm=TRUE)]
m[, n_sig_bonf := rowSums(cbind(overall_p<0.0021, knee_p<0.0021, hip_p<0.0021), na.rm=TRUE)]
# internal-replication verdict: discovery significant AND >=1 replication same-direction nominal
m[, verdict := fifelse(overall_p < 0.05 & dir_consist & (knee_p<0.05 | hip_p<0.05), "Replicated",
              fifelse(dir_consist, "Directionally consistent",
              fifelse(overall_p<0.05, "Discovery-only", "Inconsistent")))]
setorder(m, overall_p)
fwrite(m, OUTC)
cat("Wrote", OUTC, "\n"); print(m[, .(tf, overall_beta, overall_p, knee_beta, knee_p, hip_beta, hip_p, dir_consist, n_sig05, verdict)])

# ---- figure: beta across 3 OA definitions per TF, colored by significance ----
L <- d[, .(tf, outcome, est=estimate, se=se, p=p)]
L[, sig := fifelse(p < 0.0021, "Bonferroni", fifelse(p<0.05, "nominal", "ns"))]
L[, def := fifelse(outcome=="ART","overall (discovery)", fifelse(outcome=="KNEE","knee (rep.)","hip (rep.)"))]
L$def <- factor(L$def, levels=c("overall (discovery)","knee (rep.)","hip (rep.)"))
p <- ggplot(L, aes(x=tf, y=est, color=sig, group=def)) +
  geom_point(position=position_dodge(width=0.3), size=2.5) +
  geom_errorbar(aes(ymin=est-1.96*se, ymax=est+1.96*se), width=0.1, position=position_dodge(width=0.3)) +
  geom_hline(yintercept=0, linetype="dashed", color="grey40") +
  facet_wrap(~def, nrow=1) +
  scale_color_manual(values=c(Bonferroni="red", nominal="orange", ns="grey50")) +
  labs(title="Phase 2-A (interim): internal replication of TF->OA MR across OA definitions",
       subtitle="FinnGen R12 multi-SNP IVW(RE); overall=discovery, knee/hip=replications",
       x="TF", y="MR effect (genetically determined TF expression -> OA risk)") +
  theme_minimal(base_size=10) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(OUTF, p, width=9, height=4.5, dpi=150)
cat("Wrote", OUTF, "\nDONE Phase2-A(interim)\n")
