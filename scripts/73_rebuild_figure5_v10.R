# Corrected MR audit and replacement main figures for the PIEZO1-OA manuscript.
# All plotting, previewing and export are performed in R.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
  library(scales); library(svglite); library(ragg)
})
set.seed(20260831)

source("scripts/_bootstrap.R")
ROOT <- PIEZO1_ROOT
OUT  <- PIEZO1_METADATA
FIG  <- PIEZO1_FIGURES
SRC  <- file.path(PIEZO1_METADATA, "source_data")
OLD  <- reference_path("source_data")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(SRC, recursive = TRUE, showWarnings = FALSE)

pal <- c(blue="#247B9A", blue2="#73AEC4", orange="#D06A45", gold="#E2A44F",
         dark="#333333", grey="#9A9A9A", light="#E7E7E7", teal="#378C7C")
theme_pub <- function(base_size=7) theme_classic(base_size=base_size, base_family="Arial") +
  theme(axis.line=element_line(linewidth=.35), axis.ticks=element_line(linewidth=.35),
        axis.title=element_text(size=base_size), axis.text=element_text(size=base_size-.3),
        strip.background=element_blank(), strip.text=element_text(face="bold", size=base_size),
        legend.title=element_text(size=base_size-.2), legend.text=element_text(size=base_size-.4),
        plot.title=element_text(face="bold", size=base_size+.4),
        plot.tag=element_text(face="bold", size=9), panel.grid=element_blank())
theme_set(theme_pub())

save_pub <- function(p, stem, width_mm=183, height_mm=190) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite(file.path(FIG, paste0(stem,".svg")), width=w, height=h); print(p); dev.off()
  cairo_pdf(file.path(FIG, paste0(stem,".pdf")), width=w, height=h, family="Arial"); print(p); dev.off()
  agg_tiff(file.path(FIG, paste0(stem,".tiff")), width=w, height=h, units="in", res=600,
           compression="lzw"); print(p); dev.off()
  agg_png(file.path(FIG, paste0(stem,".png")), width=w, height=h, units="in", res=300);
  print(p); dev.off()
}

# ---- Correct exposure scale and MR estimators --------------------------------
snp <- fread(file.path(ROOT,"results/11_MR_snp_level.csv"))
eq <- fread(file.path(ROOT,"data/mr/eqtlgen_PIEZO1_cis.tsv"))[
  SNP %in% unique(snp$SNP), .(SNP, assessed=AssessedAllele, Z=Zscore, N=as.numeric(NrSamples))]
af <- fread(file.path(ROOT,"data/mr/eqtlgen_AF_instruments.tsv"))
af[, freqA := 1-AlleleB_all]
afl <- rbind(af[,.(SNP, assessed=AlleleA, eaf=freqA)],
             af[,.(SNP, assessed=AlleleB, eaf=AlleleB_all)])
ex <- merge(eq, afl, by=c("SNP","assessed"))
ex[, `:=`(beta_x=Z/sqrt(2*eaf*(1-eaf)*N),
          se_x=1/sqrt(2*eaf*(1-eaf)*N), F_stat=Z^2)]
snp <- merge(snp[, !c("beta_x","se_x","wald","se_wald","Fstat"), with=FALSE],
             ex[,.(SNP,beta_x,se_x,eaf_exposure=eaf,F_stat)], by="SNP")

ivw_stats <- function(r,se) {
  w <- 1/se^2; b <- sum(w*r)/sum(w); s <- sqrt(1/sum(w)); k <- length(r)
  Q <- sum(w*(r-b)^2); phi <- max(1,Q/(k-1)); sm <- s*sqrt(phi)
  C <- sum(w)-sum(w^2)/sum(w); tau2 <- max(0,(Q-k+1)/C)
  wd <- 1/(se^2+tau2); bd <- sum(wd*r)/sum(wd); sd <- sqrt(1/sum(wd))
  list(fe=c(b,s), mre=c(b,sm), dl=c(bd,sd), Q=Q, df=k-1, phi=phi, tau2=tau2)
}
wmedian <- function(r,se,B=10000) {
  w <- 1/se^2
  one <- function(x) {o<-order(x); x[o][which(cumsum(w[o])/sum(w)>=.5)[1]]}
  rb <- replicate(B, one(rnorm(length(r),r,se)))
  c(one(r), sd(rb))
}

models <- list(); weights <- list(); singles <- list(); loo <- list(); het <- list()
for (oc in unique(snp$outcome)) {
  d <- copy(snp[outcome==oc]); d[, `:=`(ratio=beta_y/beta_x, ratio_se=abs(se_y/beta_x))]
  z <- ivw_stats(d$ratio,d$ratio_se); wm <- wmedian(d$ratio,d$ratio_se)
  eg <- summary(lm(beta_y ~ beta_x, data=d, weights=1/se_y^2))
  add <- function(name,v) data.table(outcome=oc, model=name, beta=v[1], se=v[2],
    p=2*pnorm(-abs(v[1]/v[2])))
  models[[length(models)+1]] <- rbind(add("Fixed effect",z$fe),
    add("Multiplicative random effects",z$mre), add("DerSimonian-Laird",z$dl),
    add("Weighted median",wm),
    add("MR-Egger slope",c(coef(eg)[2,1],coef(eg)[2,2])))
  w <- 1/d$ratio_se^2
  weights[[length(weights)+1]] <- d[,.(outcome=oc,SNP,beta_x,se_x,F_stat,
    ratio,ratio_se,weight_pct=100*w/sum(w))]
  singles[[length(singles)+1]] <- d[,.(outcome=oc,SNP,ratio,ratio_se)]
  het[[length(het)+1]] <- data.table(outcome=oc,Q=z$Q,df=z$df,
    p=pchisq(z$Q,z$df,lower.tail=FALSE),I2=max(0,(z$Q-z$df)/z$Q)*100,phi=z$phi,tau2=z$tau2,
    egger_intercept=coef(eg)[1,1],egger_intercept_p=coef(eg)[1,4])
  for (j in seq_len(nrow(d))) {
    zz <- ivw_stats(d$ratio[-j],d$ratio_se[-j])
    loo[[length(loo)+1]] <- data.table(outcome=oc,excluded=d$SNP[j],
      beta=zz$mre[1],se=zz$mre[2],p=2*pnorm(-abs(zz$mre[1]/zz$mre[2])))
  }
}
models <- rbindlist(models); weights <- rbindlist(weights); singles <- rbindlist(singles)
loo <- rbindlist(loo); het <- rbindlist(het)
ord <- c("OA (overall, M13_ARTHROSIS)","Knee OA (gonarthrosis)","Hip OA (coxarthrosis)")
lab <- c("Overall OA","Knee OA","Hip OA")
for (x in c("models","weights","singles","loo","het")) {
  d <- get(x); d[, outcome_short:=factor(outcome,levels=ord,labels=lab)]; assign(x,d)
}
mde <- models[model %in% c("Multiplicative random effects","DerSimonian-Laird")]
mde[, `:=`(mde_beta=(qnorm(.975)+qnorm(.8))*se,
           variance_model=fifelse(model=="Multiplicative random effects","Primary","Conservative"))]
mde[, mde_OR:=exp(mde_beta)]

fwrite(ex,file.path(SRC,"58_exposure_reconstruction.csv"))
fwrite(models,file.path(SRC,"58_mr_models.csv")); fwrite(weights,file.path(SRC,"58_mr_weights.csv"))
fwrite(loo,file.path(SRC,"58_mr_leave_one_out.csv")); fwrite(het,file.path(SRC,"58_mr_diagnostics.csv"))
fwrite(mde,file.path(SRC,"58_mr_power_mde.csv"))

# ---- Figure 4 ---------------------------------------------------------------
b <- fread(file.path(OLD,"sourcedata_Fig4b_layers.csv"))
c4 <- fread(file.path(OLD,"sourcedata_Fig4c_programme_direction.csv"))
d4 <- fread(file.path(OLD,"sourcedata_Fig4d_within_state.csv"))
e4 <- fread(file.path(OLD,"sourcedata_Fig4e_kitagawa.csv"))
f4 <- fread(file.path(OLD,"sourcedata_Fig4f_piezo1_direction.csv"))

flow <- data.table(x=1:4, y=1, title=c("GSE104782","GSE57218","GSE152805","GSE51588"),
  role=c("Defines seven\nchondrocyte states","Paired cartilage\nreplication","Within-state\ncompartment support","Regional discovery;\ncross-tissue similarity"),
  col=c(pal["gold"],pal["blue"],pal["blue"],pal["grey"]))
p4a <- ggplot(flow,aes(x,y)) + geom_tile(aes(fill=col),width=.88,height=.75,colour=pal["dark"],linewidth=.35) +
  geom_text(aes(label=title),fontface="bold",vjust=-.3,size=3.1) +
  geom_text(aes(label=role),vjust=1.05,size=2.4,lineheight=.95) +
  scale_fill_identity()+coord_cartesian(xlim=c(.5,4.5),ylim=c(.48,1.52),clip="off")+theme_void()

b[, family:=factor(family,levels=c("Raw","CLR","ALR","Pairwise","GS score"))]
b[, quantity:=factor(quantity,levels=c("HomC minus preHTC contrast","HomC alone"))]
b[, ds_short:=fifelse(grepl("51588",dataset),"GSE51588 bone","GSE57218 cartilage")]
p4b <- ggplot(b,aes(family,delta,colour=ds_short)) + geom_hline(yintercept=0,colour=pal["grey"],linewidth=.3)+
  geom_point(position=position_jitterdodge(jitter.width=.12,dodge.width=.4),size=1.7,alpha=.9)+
  facet_wrap(~quantity,ncol=1,scales="free_y")+
  scale_colour_manual(values=c("GSE51588 bone"=unname(pal["blue"]),"GSE57218 cartilage"=unname(pal["orange"])))+
  labs(x=NULL,y="Layer-specific estimate",colour=NULL)+theme(legend.position="top",strip.placement="outside")

c4[, programme:=factor(programme,levels=c("EC","HomC","HTC","RegC","FC","preHTC","ProC"))]
p4c <- ggplot(c4,aes(factor(donor),programme,fill=dir)) + geom_tile(colour="white",linewidth=.8)+
  geom_text(aes(label=ifelse(dir=="lateral higher","+","-")),colour="white",fontface="bold",size=3.5)+
  scale_fill_manual(values=c("lateral higher"=unname(pal["blue"]),"medial higher"=unname(pal["orange"])))+
  labs(x="Donor",y=NULL,fill=NULL,title="Programme direction")+theme(legend.position="none")

d4[, state:=factor(state,levels=c("EC","HomC","HTC","RegC","FC","preHTC","ProC"))]
p4d <- ggplot(d4,aes(factor(donor),state,fill=delta_lateral_minus_medial))+
  geom_tile(colour="white",linewidth=.8)+geom_text(aes(label=sprintf("%+.02f",delta_lateral_minus_medial)),size=2.2)+
  scale_fill_gradient2(low=pal["orange"],mid="white",high=pal["blue"],midpoint=0)+
  labs(x="Donor",y=NULL,fill="Lateral - medial",title="Within-state PIEZO1")+theme(legend.position="none")

el <- melt(e4,id.vars=c("donor","total_delta_lateral_minus_medial","pct_within"),
  measure.vars=c("within_state_component","composition_component"),variable.name="component")
p4e <- ggplot(el,aes(factor(donor),value,fill=component)) + geom_col()+
  geom_point(data=e4,aes(factor(donor),total_delta_lateral_minus_medial),inherit.aes=FALSE,shape=18,size=2.2)+
  scale_fill_manual(values=c("within_state_component"=unname(pal["blue"]),"composition_component"=unname(pal["light"])),
    labels=c("Within state","Composition"))+geom_hline(yintercept=0,linewidth=.3)+
  labs(x="Donor",y="Contribution",fill=NULL,title="Decomposition")+theme(legend.position="none")

f4[, ylab_clean:=fcase(src=="GSE104782","GSE104782: correlation with state axis\n10 donors",
                       src=="GSE152805","GSE152805: lateral - medial mean\n3 donors",
                       src=="GSE57218","GSE57218: preserved - OA-affected\n33 donors",
                       default="GSE51588: lateral - medial, OA\n20 donors")]
f4[, ylab_clean:=factor(ylab_clean,levels=rev(unique(ylab_clean)))]
p4f <- ggplot(f4,aes(est,ylab_clean,colour=tissue))+geom_vline(xintercept=0,linetype=2,colour=pal["grey"])+
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=.13,na.rm=TRUE)+geom_point(size=2.2)+
  scale_colour_manual(values=c("Cartilage"=unname(pal["orange"]),"Subchondral bone"=unname(pal["blue"])))+
  labs(x="Dataset-specific estimate (95% CI where available)",y=NULL,colour=NULL,
       title="PIEZO1 direction is context dependent")+theme(legend.position="top")

fig4 <- p4a / p4b / (p4c|p4d|p4e) / p4f + plot_layout(heights=c(.55,1.35,1.25,1.05)) +
  plot_annotation(tag_levels="a")
save_pub(fig4,"Fig4_programme_contrast_PIEZO1_direction_corrected",183,205)

# ---- Figure 5 ---------------------------------------------------------------
models[, model:=factor(model,levels=c("Fixed effect","Multiplicative random effects","DerSimonian-Laird",
  "Weighted median","MR-Egger slope"))]
forest <- function(dd, title, highlight_primary = FALSE) {
  p <- ggplot(dd, aes(exp(beta), model, colour = model))
  if (highlight_primary) {
    p <- p + annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1.5, ymax = 2.5,
                      fill = pal["blue2"], alpha = 0.13,
                      colour = pal["blue"], linewidth = 0.55)
  }
  p + geom_vline(xintercept = 1, linetype = 2, colour = pal["grey"]) +
    geom_errorbarh(aes(xmin = exp(beta - 1.96 * se), xmax = exp(beta + 1.96 * se)), height = .12) +
    geom_point(size = 2) + facet_wrap(~outcome_short, nrow = 1) +
    scale_colour_manual(values = unname(c(pal["grey"], pal["blue"], pal["orange"], pal["teal"], pal["gold"]))) +
    labs(x = "Odds ratio per 1-SD higher predicted blood PIEZO1", y = NULL,
         title = title,
         subtitle = if (highlight_primary) "Outlined row: primary multiplicative random-effects model" else NULL) +
    theme(legend.position = "none")
}
p5a <- forest(models[model %in% levels(model)[1:3]],
              "IVW estimates across variance models", highlight_primary = TRUE)
p5b <- forest(models[model %in% levels(model)[4:5]], "Sensitivity estimators")

wplot <- merge(singles,weights[,.(outcome,SNP,weight_pct)],by=c("outcome","SNP"))
wplot[, SNP:=factor(SNP,levels=rev(unique(SNP[order(weight_pct)])))]
p5c <- ggplot(wplot,aes(exp(ratio),SNP,size=weight_pct))+
  geom_vline(xintercept=1,linetype=2,colour=pal["grey"])+
  geom_errorbarh(aes(xmin=exp(ratio-1.96*ratio_se),xmax=exp(ratio+1.96*ratio_se)),height=.1,colour=pal["grey"])+
  geom_point(colour=pal["blue"])+facet_wrap(~outcome_short,nrow=1,scales="free_x")+
  scale_size_area(max_size=6)+labs(x="Per-SNP Wald ratio",y=NULL,size="IVW weight (%)",
    title="Single-instrument estimates")+theme(legend.position="top")

co <- fread(file.path(ROOT,"results/41_P0_coloc_PPH3_report.csv"))
co[, outcome_short:=factor(outcome,levels=c("OA (overall)","Knee OA","Hip OA"),labels=lab)]
cl <- melt(co,id.vars=c("outcome_short","n_regsig_fg","lead_distance_kb"),
  measure.vars=c("PP.H3","PP.H4"),variable.name="hypothesis",value.name="posterior")
p5d <- ggplot(cl,aes(outcome_short,posterior,colour=hypothesis,group=hypothesis))+
  geom_line(linewidth=.5)+geom_point(size=2.2)+scale_y_log10(labels=label_percent(accuracy=.001),limits=c(1e-5,1.2))+
  scale_colour_manual(values=c("PP.H3"=unname(pal["blue"]),"PP.H4"=unname(pal["gold"])),labels=c("Distinct signals (H3)","Shared variant (H4)"))+
  labs(x=NULL,y="Posterior probability (log scale)",colour=NULL,title="Colocalization by outcome")+
  theme(legend.position="top",axis.text.x=element_text(angle=20,hjust=1))

p5e <- ggplot(mde,aes(outcome_short,mde_OR,fill=variance_model))+
  geom_hline(yintercept=1,linewidth=.3)+geom_col(position=position_dodge(.7),width=.62)+
  geom_text(aes(label=sprintf("%.2f",mde_OR)),position=position_dodge(.7),vjust=-.4,size=2.5)+
  scale_fill_manual(values=c("Primary"=unname(pal["blue"]),"Conservative"=unname(pal["orange"])))+
  coord_cartesian(ylim=c(1,1.17))+labs(x=NULL,y="Minimum detectable OR (80% power)",fill=NULL,title="Power boundary")+
  theme(legend.position="top",axis.text.x=element_text(angle=20,hjust=1))

fig5 <- p5a / p5b / p5c / (p5d|p5e) + plot_layout(heights=c(1, .9, 1.2, 1)) +
  plot_annotation(tag_levels="a")
save_pub(fig5,"Figure5_MR_coloc_power_v11",183,195)

# Corrected supplementary audit: the valid F statistic is Z^2; old MAF-based F double-counted genotype variance.
ivp <- ex[,.(SNP,F_stat,eaf,beta_x,se_x)]
p13a <- ggplot(ivp,aes(reorder(SNP,F_stat),F_stat))+geom_col(fill=pal["blue"],width=.65)+
  geom_hline(yintercept=10,linetype=2)+coord_flip()+scale_y_log10()+
  labs(x=NULL,y="F statistic (log scale)",title="Instrument strength: F = Z^2")
p13b <- ggplot(weights[outcome==ord[1]],aes(reorder(SNP,weight_pct),weight_pct))+geom_col(fill=pal["orange"],width=.65)+
  geom_text(aes(label=sprintf("%.1f%%",weight_pct)),hjust=-.1,size=2.5)+coord_flip(clip="off")+
  scale_y_continuous(limits=c(0,95),expand=expansion(mult=c(0,.03)))+
  labs(x=NULL,y="IVW weight (%)",title="Overall-OA weight concentration")
p13c <- ggplot(loo,aes(exp(beta),excluded,colour=outcome_short))+
  geom_vline(xintercept=1,linetype=2)+geom_errorbarh(aes(xmin=exp(beta-1.96*se),xmax=exp(beta+1.96*se)),height=.1)+geom_point()+
  facet_wrap(~outcome_short,ncol=1,scales="free_x")+labs(x="Leave-one-out OR",y="Excluded SNP",colour=NULL,title="Leave-one-out sensitivity")+
  theme(legend.position="none")
figS13 <- (p13a|p13b|p13c)+plot_layout(widths=c(1,1,1.25))+plot_annotation(tag_levels="a")
save_pub(figS13,"FigS13_corrected_MR_scale_audit",183,85)

# Recompute external validation using allele frequencies for the validation variants.
gv <- fread(file.path(ROOT,"results/50_gtex_scale_validation.csv"))
afall <- fread(file.path(ROOT,"data/mr/eqtlgen_AF.txt.gz"),select=c("SNP","AlleleA","AlleleB","AlleleB_all"))
afall <- afall[SNP %in% gv$rsid]; afall[,freqA:=1-AlleleB_all]
gaf <- rbind(afall[,.(rsid=SNP,allele=AlleleA,eaf=freqA)],afall[,.(rsid=SNP,allele=AlleleB,eaf=AlleleB_all)])
gv <- merge(gv,gaf,by.x=c("rsid","eQTLGen_assessed"),by.y=c("rsid","allele"))
gv[, beta_corrected:=Z/sqrt(2*eaf*(1-eaf)*N_eQTLGen)]
fit <- lm(GTEx_nes_aligned~beta_corrected,data=gv)
p14 <- ggplot(gv,aes(beta_corrected,GTEx_nes_aligned,colour=abs(Z)>=15))+geom_hline(yintercept=0,linetype=2,colour=pal["grey"])+
  geom_vline(xintercept=0,linetype=2,colour=pal["grey"])+geom_abline(slope=1,linetype=3)+geom_smooth(method="lm",se=FALSE,colour=pal["dark"],linewidth=.5)+
  geom_point(size=1.8)+scale_colour_manual(values=unname(c(pal["blue2"],pal["orange"])),labels=c("|Z| < 15","|Z| >= 15"))+
  annotate("text",x=-Inf,y=Inf,hjust=-.05,vjust=1.2,label=sprintf("slope = %.2f; r = %.2f; n = %d",coef(fit)[2],cor(gv$beta_corrected,gv$GTEx_nes_aligned),nrow(gv)),size=2.8)+
  labs(x="eQTLGen per-allele effect reconstructed with allele frequency",y="GTEx whole-blood NES",colour=NULL)+theme(legend.position="top")
save_pub(p14,"FigS14_corrected_exposure_scale_validation",100,85)
fwrite(gv,file.path(SRC,"58_gtex_corrected_scale_validation.csv"))

cat("Completed corrected MR and figure exports.\n")
