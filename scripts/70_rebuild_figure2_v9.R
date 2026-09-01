suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(svglite); library(ragg)
})

source("scripts/_bootstrap.R")
SRC <- reference_path("source_data")
FIG <- PIEZO1_FIGURES
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

pal <- c(blue="#2A819E", blue_light="#C8DDE5", orange="#D46A45", gold="#D8A43A",
         grey="#9B9B9B", light="#F2F2F2", dark="#2B2B2B")
theme_pub <- function(base_size=7.5) theme_classic(base_size=base_size, base_family="Arial") +
  theme(axis.line=element_line(linewidth=.35), axis.ticks=element_line(linewidth=.35),
        axis.title=element_text(size=base_size), axis.text=element_text(size=base_size-.2),
        plot.title=element_text(size=base_size+.5,face="bold",margin=margin(b=4)),
        plot.subtitle=element_text(size=base_size-.3,colour="#666666",margin=margin(b=5)),
        plot.tag=element_text(size=9,face="bold"), strip.background=element_blank(),
        strip.text=element_text(size=base_size,face="bold"), panel.grid=element_blank())
theme_set(theme_pub())

save_pub <- function(p, stem, width_mm, height_mm) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite(file.path(FIG,paste0(stem,".svg")),width=w,height=h); print(p); dev.off()
  cairo_pdf(file.path(FIG,paste0(stem,".pdf")),width=w,height=h,family="Arial"); print(p); dev.off()
  agg_tiff(file.path(FIG,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw"); print(p); dev.off()
  agg_png(file.path(FIG,paste0(stem,".png")),width=w,height=h,units="in",res=300); print(p); dev.off()
}

# Figure 2: paired compartment discovery, transcriptome context, pathways, and model hierarchy.
pa <- fread(file.path(SRC,"sourcedata_Fig2a_paired_piezo1.csv"))
setnames(pa,c("donor_id","compartment","expr"),c("donor","region","expression"))
pa[, region := factor(region, levels=c("Lateral","Medial"))]
p2a <- ggplot(pa,aes(region,expression,group=donor)) +
  geom_line(colour="#C7C7C7",linewidth=.35) +
  geom_point(aes(fill=region),shape=21,size=1.9,colour=pal["dark"],stroke=.25) +
  stat_summary(aes(group=1),fun=mean,geom="line",linewidth=.9,colour=pal["dark"]) +
  stat_summary(aes(group=1),fun=mean,geom="point",shape=18,size=3,colour=pal["dark"]) +
  scale_fill_manual(values=c(Lateral="#2A819E",Medial="#D46A45"),guide="none") +
  labs(x=NULL,y="PIEZO1 expression",title="Paired compartment contrast",
       subtitle="Medial - lateral = -0.504; limma P=0.000273")

vol <- fread(file.path(SRC,"Figure2B_paired_volcano_source.csv"))
vol[, sig := adj.P.Val < .05]
vol[, y := -log10(pmax(P.Value,1e-300))]
p2b <- ggplot(vol,aes(logFC,y)) +
  geom_point(aes(colour=sig),size=.45,alpha=.62) +
  geom_vline(xintercept=0,colour="#A0A0A0",linewidth=.3) +
  geom_point(data=vol[is_piezo==TRUE],shape=21,size=2.7,fill=pal["gold"],colour=pal["dark"],stroke=.45) +
  geom_text(data=vol[is_piezo==TRUE],aes(label="PIEZO1"),nudge_x=-.15,nudge_y=.8,hjust=1,size=2.5) +
  scale_colour_manual(values=c(`FALSE`="#D7DADD",`TRUE`="#D46A45"),guide="none") +
  labs(x="Medial minus lateral log2 fold change",y="-log10(P)",
       title="Paired transcriptome contrast",subtitle="PIEZO1 highlighted; orange indicates FDR < 0.05")

go <- fread(file.path(SRC,"sourcedata_Fig2c_go.csv"))
if (!"side" %in% names(go)) go[, side := ifelse(grepl("medial",direction,ignore.case=TRUE),"Medial","Lateral")]
if (!"score" %in% names(go)) go[, score := ifelse(side=="Medial",1,-1) * -log10(p.adjust)]
go[, Description := factor(Description,levels=Description[order(score)])]
p2c <- ggplot(go,aes(score,Description)) +
  geom_vline(xintercept=0,colour="#A0A0A0",linewidth=.35) +
  geom_segment(aes(x=0,xend=score,y=Description,yend=Description,colour=side),linewidth=.7) +
  geom_point(aes(size=Count,fill=side),shape=21,colour="white",stroke=.35) +
  scale_colour_manual(values=c(Lateral="#2A819E",Medial="#D46A45"),guide="none") +
  scale_fill_manual(values=c(Lateral="#2A819E",Medial="#D46A45"),guide="none") +
  scale_size(range=c(2.2,4.5),name="Genes") +
  labs(x="Signed -log10(FDR)",y=NULL,title="Opposing regional programmes",
       subtitle="Blue: lateral enriched; orange: medial enriched") +
  theme(legend.position="bottom",legend.key.height=grid::unit(3,"mm"))

md <- fread(file.path(SRC,"sourcedata_Fig2d_models.csv"))
md[lab == "OA donors (n = 20)", P := 0.000272626386996261]
md[, lab := factor(lab,levels=rev(lab))]
md[, p_label := paste0("P=",fcase(P<.01,formatC(P,format="f",digits=5),
                                  P<.1,formatC(P,format="f",digits=3),
                                  default=formatC(P,format="f",digits=2)))]
p2d <- ggplot(md,aes(estimate,lab)) +
  geom_vline(xintercept=0,linetype=2,colour="#8E8E8E",linewidth=.4) +
  geom_errorbarh(aes(xmin=ci_lo,xmax=ci_hi,colour=emph),height=.12,linewidth=.6) +
  geom_point(aes(fill=emph),shape=21,size=2.5,colour=pal["dark"],stroke=.35) +
  geom_label(aes(x=ci_hi+.07,label=p_label),hjust=0,size=2.3,colour="#555555",
             fill="white",label.size=NA,label.padding=grid::unit(.06,"lines")) +
  scale_colour_manual(values=c(primary="#D46A45",other="#9B9B9B"),guide="none") +
  scale_fill_manual(values=c(primary="#D46A45",other="#D6D6D6"),guide="none") +
  labs(x="Estimate (95% CI)",y=NULL,title="Model hierarchy and sensitivity") +
  coord_cartesian(xlim=c(min(md$ci_lo)-.08,max(md$ci_hi)+.68),clip="off") +
  theme(axis.text.y=element_text(size=6.5,lineheight=.92))

fig2 <- (p2a|p2b) / (p2c|p2d) +
  plot_layout(heights=c(1,1.05),widths=c(.92,1.18)) + plot_annotation(tag_levels="a")
save_pub(fig2,"Figure2_subchondral_compartment_v9",183,145)
fwrite(md, file.path(FIG, "SourceData_Figure2d_models_v9.csv"))

# Supplementary Figure S1: move annotations into dedicated right-hand text space.
ds <- fread(file.path(SRC,"Figure1_disease_level_effects.csv"))
ds[, label := factor(label,levels=rev(label))]
pS1a <- ggplot(ds,aes(Diff,label)) +
  geom_vline(xintercept=0,linetype=2,colour="#999999",linewidth=.4) +
  geom_errorbarh(aes(xmin=CI_low,xmax=CI_high,colour=Tissue),height=.12,linewidth=.65) +
  geom_point(aes(fill=Tissue),shape=21,size=2.5,colour="white",stroke=.35) +
  geom_text(aes(x=CI_high+.12,label=paste0("P=",format(P,digits=2))),hjust=0,size=2.45,colour="#555555") +
  scale_colour_manual(values=c("Subchondral bone"="#D8A43A","Synovium"="#4D9F94"),guide="none") +
  scale_fill_manual(values=c("Subchondral bone"="#D8A43A","Synovium"="#4D9F94"),guide="none") +
  labs(x="OA - control expression difference (95% CI)",y=NULL,title="No uniform disease-level direction") +
  coord_cartesian(xlim=c(min(ds$CI_low)-.12,max(ds$CI_high)+.72),clip="off")

au <- fread(file.path(SRC,"SupplementaryFigureS1_auc.csv"))
au[, label := factor(label,levels=rev(label))]
pS1b <- ggplot(au,aes(AUC,label)) +
  geom_vline(xintercept=.5,linetype=2,colour="#999999",linewidth=.4) +
  geom_segment(aes(x=.5,xend=AUC,y=label,yend=label),linewidth=.8,colour="#AAB2B7") +
  geom_point(aes(fill=AUC>.8),shape=21,size=2.6,colour="white",stroke=.35) +
  geom_text(aes(x=AUC+.018,label=sprintf("%.2f",AUC)),hjust=0,size=2.45,colour="#555555") +
  scale_fill_manual(values=c(`FALSE`="#9EABB2",`TRUE`=pal["gold"]),guide="none") +
  labs(x="Single-gene AUC",y=NULL,title="Discrimination is inconsistent") +
  coord_cartesian(xlim=c(.48,.88),clip="off")

figS1 <- (pS1a | pS1b) + plot_layout(widths=c(1.15,1)) + plot_annotation(tag_levels="a")
save_pub(figS1,"SuppFigS1_disease_context_fixed",183,78)

cat("Fixed Figure 2 and Supplementary Figure S1 written to",FIG,"\n")
