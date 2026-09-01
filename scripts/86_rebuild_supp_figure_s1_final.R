suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

source("scripts/_bootstrap.R")
SRC <- reference_path("source_data")
FIG <- PIEZO1_FIGURES
pal <- c(gold="#D8A43A", dark="#2B2B2B")
theme_set(theme_classic(base_size=7.5, base_family="Arial") +
  theme(axis.line=element_line(linewidth=.35), axis.ticks=element_line(linewidth=.35),
        axis.title=element_text(size=7.5), axis.text=element_text(size=7.3),
        plot.title=element_text(size=8,face="bold",margin=margin(b=4)),
        plot.tag=element_text(size=9,face="bold"), panel.grid=element_blank()))
save_pub <- function(p, stem, width_mm, height_mm) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite(file.path(FIG,paste0(stem,".svg")),width=w,height=h); print(p); dev.off()
  cairo_pdf(file.path(FIG,paste0(stem,".pdf")),width=w,height=h,family="Arial"); print(p); dev.off()
  agg_tiff(file.path(FIG,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw"); print(p); dev.off()
  agg_png(file.path(FIG,paste0(stem,".png")),width=w,height=h,units="in",res=300); print(p); dev.off()
}

ds <- fread(file.path(SRC,"Figure1_disease_level_effects.csv"))
ds[, label := factor(label,levels=rev(label))]
ds[, p_label := ifelse(grepl("GSE51588", as.character(label)), "P=0.500", sprintf("P=%.3f", P))]
pS1a <- ggplot(ds,aes(Diff,label)) +
  geom_vline(xintercept=0,linetype=2,colour="#999999",linewidth=.4) +
  geom_errorbarh(aes(xmin=CI_low,xmax=CI_high,colour=Tissue),height=.12,linewidth=.65) +
  geom_point(aes(fill=Tissue),shape=21,size=2.5,colour="white",stroke=.35) +
  geom_text(aes(x=CI_high+.12,label=p_label),hjust=0,size=2.45,colour="#555555") +
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
save_pub(figS1,"Figure_S1_Disease_Context_P0500",183,78)
