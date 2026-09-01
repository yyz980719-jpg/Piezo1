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
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
pal <- c(blue="#2A819E", blue_light="#C8DDE5", orange="#D46A45", gold="#E6AA50",
         grey="#9B9B9B", light="#F2F2F2", dark="#2B2B2B", white="#FFFFFF")
theme_pub <- function(base_size=7) theme_classic(base_size=base_size, base_family="Arial") +
  theme(axis.line=element_line(linewidth=.35), axis.ticks=element_line(linewidth=.35),
        axis.title=element_text(size=base_size), axis.text=element_text(size=base_size-.2),
        plot.title=element_text(size=base_size+.5,face="bold"),
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
wrap_txt <- function(x,width) vapply(x,function(z) paste(strwrap(z,width=width),collapse="\n"),character(1))

arch <- data.table(
  y=6:1,
  step=c("Disease-level expression","Paired joint compartments","Donor-aware states",
         "Cross-dataset concordance","Within-state decomposition","Blood-expression genetics"),
  data=c("Four tissue datasets","GSE51588: 20 OA + 5 non-OA donors",
         "GSE104782: 1,464 cells; 10 donors","GSE57218: 33 pairs; GSE152805: 3 pairs",
         "GSE152805: 7 states x 3 donors","eQTLGen cis-eQTL x FinnGen R12"),
  unit=c("Sample within dataset","Paired donor difference","One estimate per donor",
         "Paired donor difference","Donor x state mean","Instrument-level Wald ratio"),
  family=c("context","regional","regional","regional","decomposition","genetics"))
arch[, fill:=fcase(family=="context",pal["grey"],family=="decomposition",pal["gold"],
                   family=="genetics",pal["orange"],default=pal["blue"])]
p1a <- ggplot(arch) +
  geom_rect(aes(xmin=0,xmax=1,ymin=y-.38,ymax=y+.38),fill=pal["light"],colour=pal["dark"],linewidth=.35)+
  geom_rect(aes(xmin=0,xmax=.018,ymin=y-.38,ymax=y+.38,fill=fill),colour=NA)+scale_fill_identity()+
  geom_text(aes(x=.045,y=y,label=wrap_txt(step,22)),hjust=0,fontface="bold",size=2.7,lineheight=.95)+
  geom_text(aes(x=.31,y=y,label=wrap_txt(data,31)),hjust=0,colour=pal["dark"],size=2.45,lineheight=.95)+
  geom_text(aes(x=.71,y=y,label=wrap_txt(unit,28)),hjust=0,colour="#666666",fontface="italic",size=2.35,lineheight=.95)+
  annotate("text",x=c(.045,.31,.71),y=6.67,label=c("Evidence layer","Data and sample structure","Inferential unit"),
           hjust=0,fontface="bold",size=2.65)+
  coord_cartesian(xlim=c(0,1),ylim=c(.5,6.9),clip="off")+theme_void(base_family="Arial")

d1 <- fread(file.path(SRC,"sourcedata_Fig1b_disease_contrasts.csv"))
d1[, label:=factor(paste0(Dataset,"\n",Tissue),levels=rev(paste0(Dataset,"\n",Tissue)))]
d1[, p_label := ifelse(Dataset == "GSE51588", "P=0.500", sprintf("P=%.3f", P))]
p1b <- ggplot(d1,aes(Diff,label)) + geom_vline(xintercept=0,linetype=2,colour=pal["grey"])+
  geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.12,colour=pal["grey"],linewidth=.55)+
  geom_point(size=2.5,fill=pal["blue_light"],colour=pal["dark"],shape=21)+
  geom_text(aes(x=CI_high+.06,label=p_label),hjust=0,size=2.3,colour="#666666")+
  labs(x="PIEZO1 difference: OA minus control (95% CI)",y=NULL,title="Disease-level contrasts are inconsistent")+
  coord_cartesian(xlim=c(min(d1$CI_low)-.1,max(d1$CI_high)+.45),clip="off")

ct <- fread(file.path(SRC,"sourcedata_Fig1c_interpretation_contract.csv"))
ct[layer == "Paired joint compartments", can := "Compartment-associated expression within OA joints, free of between-person confounding"]
ct[layer == "Donor-aware chondrocyte states", cannot := "Cell-level independence, or evidence that the state transition is causal"]
ct[layer == "Programme-score replication", `:=`(
  layer = "Programme-score concordance",
  can = "A consistent programme-level contrast across the analysed paired datasets"
)]
ct[layer == "Within-state decomposition", cannot := "Three donors provide directional support only, not independent validation"]
ct[layer == "Blood-expression Mendelian randomization", `:=`(
  can = "Bounds for a systemic effect under the whole-blood instrument model",
  cannot = "Any local joint-tissue mechanism, or tissue-specific expression component"
)]
ct[, y:=rev(seq_len(.N))]
ct[, `:=`(layer=wrap_txt(layer,16), supports=wrap_txt(can,21), limits=wrap_txt(cannot,21))]
p1c <- ggplot(ct) +
  geom_rect(aes(xmin=0,xmax=1,ymin=y-.43,ymax=y+.43,fill=factor(y%%2)),colour="#D0D0D0",linewidth=.25)+
  scale_fill_manual(values=c("0"="#FFFFFF","1"="#F7F7F7"),guide="none")+
  geom_vline(xintercept=c(.27,.635),colour="#C8C8C8",linewidth=.3)+
  geom_text(aes(x=.015,y=y,label=layer),hjust=0,fontface="bold",size=2.05,lineheight=.9)+
  geom_text(aes(x=.285,y=y,label=supports),hjust=0,colour=pal["blue"],size=1.92,lineheight=.9)+
  geom_text(aes(x=.65,y=y,label=limits),hjust=0,colour=pal["orange"],size=1.92,lineheight=.9)+
  annotate("text",x=c(.015,.285,.65),y=max(ct$y)+.68,label=c("Evidence layer","SUPPORTS","DOES NOT SUPPORT"),
           hjust=0,fontface="bold",colour=c(pal["dark"],pal["blue"],pal["orange"]),size=2.55)+
  coord_cartesian(xlim=c(0,1),ylim=c(.45,max(ct$y)+.85),clip="off")+theme_void(base_family="Arial")

fig1 <- p1a / (p1b|p1c) + plot_layout(heights=c(1.02,1.18),widths=c(.78,1.52)) +
  plot_annotation(tag_levels="a")
save_pub(fig1,"Figure1_Final_P0500",183,170)
