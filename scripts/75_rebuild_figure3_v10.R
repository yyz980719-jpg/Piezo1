suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(svglite); library(ragg)
})

source("scripts/_bootstrap.R")
OUT <- PIEZO1_METADATA
SRC <- reference_path("source_data")
FIG <- PIEZO1_FIGURES
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)

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

# Figure 1 contract: separate evidence roles from claims that the data cannot support.
arch <- data.table(
  y=6:1,
  step=c("Disease-level expression","Paired joint compartments","Donor-aware states",
         "Cross-dataset replication","Within-state decomposition","Blood-expression genetics"),
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
d1[, label:=paste0(Dataset,"\n",Tissue)]
d1[, label:=factor(label,levels=rev(label))]
p1b <- ggplot(d1,aes(Diff,label)) + geom_vline(xintercept=0,linetype=2,colour=pal["grey"])+
  geom_errorbarh(aes(xmin=CI_low,xmax=CI_high),height=.12,colour=pal["grey"],linewidth=.55)+
  geom_point(size=2.5,fill=pal["blue_light"],colour=pal["dark"],shape=21)+
  geom_text(aes(x=CI_high+.06,label=paste0("P=",format(P,digits=2))),hjust=0,size=2.3,colour="#666666")+
  labs(x="PIEZO1 difference: OA minus control (95% CI)",y=NULL,title="Disease-level contrasts are inconsistent")+
  coord_cartesian(xlim=c(min(d1$CI_low)-.1,max(d1$CI_high)+.45),clip="off")

ct <- fread(file.path(SRC,"sourcedata_Fig1c_interpretation_contract.csv"))
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
save_pub(fig1,"Fig1_evidence_architecture_fixed",183,170)

# Figure 3 contract: donor-level association is primary; categorical state contrasts remain secondary.
axis_df <- data.table(x=1:7,state=c("ProC","EC","HomC","preHTC","HTC","RegC","FC"))
p3a <- ggplot(axis_df,aes(x,1))+
  geom_segment(aes(x=1,xend=7,y=1,yend=1),arrow=arrow(length=grid::unit(2,"mm")),linewidth=.7,colour=pal["dark"])+
  geom_point(aes(fill=x),shape=21,size=4,colour="white",stroke=.5)+
  geom_text(aes(y=.83,label=state),size=2.5,fontface="bold")+
  annotate("text",x=1,y=1.18,label="homeostatic pole",hjust=0,colour=pal["blue"],fontface="bold",size=2.4)+
  annotate("text",x=7,y=1.18,label="hypertrophic pole",hjust=1,colour=pal["orange"],fontface="bold",size=2.4)+
  annotate("text",x=4,y=.62,label="Ordered by median marker-axis score; schematic, not a mechanistic trajectory",size=2.25,colour="#666666")+
  scale_fill_gradient(low=pal["blue"],high=pal["orange"],guide="none")+
  coord_cartesian(xlim=c(.7,7.3),ylim=c(.5,1.28),clip="off")+theme_void(base_family="Arial")

rho <- fread(file.path(SRC,"sourcedata_Fig3b_donor_rho.csv"))
axis_q <- fread(file.path(OUT, "03_AxisDonor_Q1Q4.csv"))
rho[, c("mean_Q1_homeo", "mean_Q4_hyp") := NULL]
rho <- merge(rho, axis_q[, .(patient,
                             mean_Q1_homeo = Q1,
                             mean_Q4_hyp = Q4)], by = "patient")
rho[, patient:=factor(patient,levels=patient[order(rho)] )]
p3b <- ggplot(rho,aes(rho,patient)) + geom_vline(xintercept=0,linewidth=.4)+
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=.1,colour="#C8C8C8",linewidth=.55)+
  geom_point(aes(fill=rho>0),shape=21,size=2.5,colour=pal["dark"])+
  scale_fill_manual(values=c("FALSE"="#2A819E","TRUE"="#D46A45"),guide="none")+
  annotate("rect",xmin=-.216,xmax=-.058,ymin=.55,ymax=1.05,fill=pal["blue_light"],colour=NA)+
  annotate("point",x=-.137,y=.8,shape=18,size=3,colour=pal["blue"])+
  labs(x="Spearman rho: PIEZO1 versus state axis",y="Donor",title="Nine of ten donor estimates are negative")+
  coord_cartesian(xlim=c(min(rho$lo)-.02,max(rho$hi)+.03),clip="off")

pair <- melt(rho[,.(patient,Q1=mean_Q1_homeo,Q4=mean_Q4_hyp)],id.vars="patient",variable.name="quartile",value.name="value")
pair[, quartile:=factor(quartile,levels=c("Q1","Q4"),labels=c("Q1 most\nhomeostatic","Q4 most\nhypertrophic"))]
p3c <- ggplot(pair,aes(quartile,value,group=patient))+
  geom_line(colour="#C7C7C7",linewidth=.45)+geom_point(aes(fill=quartile),shape=21,size=2.2,colour=pal["dark"])+
  stat_summary(aes(group=1),fun=mean,geom="line",linewidth=1,colour=pal["dark"])+
  stat_summary(aes(group=1),fun=mean,geom="point",shape=18,size=3,colour=pal["dark"])+
  scale_fill_manual(values=c("Q1 most\nhomeostatic"="#2A819E","Q4 most\nhypertrophic"="#D46A45"),guide="none")+
  labs(x=NULL,y="Mean PIEZO1 (donor level)",title="Donor-within quartiles: paired t P=0.013")

pb <- fread(file.path(SRC,"sourcedata_Fig3d_pseudobulk_contrasts.csv"))
pb[, label:=gsub("_minus_"," - ",contrast)]
pb[, label:=factor(label,levels=rev(label))]
p3d <- ggplot(pb,aes(logFC,label,colour=sig))+
  geom_vline(xintercept=0,linewidth=.4)+geom_errorbarh(aes(xmin=CI.L,xmax=CI.R),height=.12,linewidth=.55)+geom_point(size=2.5)+
  geom_text(aes(x=CI.R+.08,label=paste0("FDR=",format(FDR,digits=2))),hjust=0,size=2.2,colour="#666666")+
  scale_colour_manual(values=c("FDR < 0.05"="#D46A45","not significant"="#79AFC0"),guide="none")+
  labs(x="PIEZO1 log2 fold change (95% CI)",y=NULL,title="Donor-blocked pseudobulk contrasts")+
  coord_cartesian(xlim=c(min(pb$CI.L)-.1,max(pb$CI.R)+.55),clip="off")

hm <- fread(file.path(SRC,"sourcedata_Fig3e_heatmap.csv"))
hm[, patient:=factor(patient,levels=rev(c("OA1","OA10","OA2","OA3","OA4","OA5","OA6","OA7","OA8","OA9")))]
hm[, state:=factor(state,levels=c("HTC","FC","preHTC","EC","RegC","HomC","ProC"))]
p3e <- ggplot(hm,aes(state,patient,fill=z))+geom_tile(colour="white",linewidth=.5)+
  scale_fill_gradient2(low=pal["blue"],mid="#F4F4F4",high=pal["orange"],midpoint=0,limits=c(-2,2),na.value="#D0D0D0")+
  labs(x=NULL,y="Donor",fill="Within-donor\nz score",title="State-level expression pattern")+
  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="right")

fig3 <- (p3a|p3b) / (p3c|p3d|p3e) + plot_layout(heights=c(.72,1.15),widths=c(1,1,1)) +
  plot_annotation(tag_levels="a")
save_pub(fig3,"Figure3_donor_aware_state_axis_v10",183,135)

cat("Fixed Figure 1 and Figure 3 written to",FIG,"\n")
