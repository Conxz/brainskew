library(ggpubr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(ggrepel)
library(tidyr)
library(GGally)
library(plyr)
library(xtable)
options(stringsAsFactors = FALSE)
#----------------------------------------------------------------------
# define working_dir
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
working_dir=paste(dir,"lg-ukbiobank/analysis/xiangzhen/brainTorque/torque_ukb_N40681/GWAS_Dick/calc_h2/reml/",sep="")
out_dir=paste(dir,"lg-ukbiobank/analysis/xiangzhen/brainTorque/torque_ukb_N40681/GWAS_Dick/calc_h2/output/",sep="")
dir.create(file.path(out_dir), showWarnings = FALSE)
setwd(working_dir)

# define pattern for this run
#pattern_run2="_noBioCovs_noAssessmentC"
#pattern_run2="_noBioCovs_noAssessmentC_TBV"
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Univariate REML runs
tables<-list.files(working_dir,pattern=".table")

for (t in tables){
  tmp<-read.table(t)
  if (t==tables[1]) { summary_h2 <-tmp} else { summary_h2 <-rbind(summary_h2,tmp)}  
  rm(tmp)
};rm(t)
rm(tables)
# extract relevant info; organize into cols
colnames(summary_h2)<-c("file","n","h2","se","pval")
summary_h2$file<-gsub(":n","",summary_h2$file)
summary_h2$phenotype<-sapply(strsplit(summary_h2$file,"_"),"[[",1)
summary_h2$Dist<-sapply(strsplit(summary_h2$file,"_"),"[[",5)
summary_h2$Asy<-sapply(strsplit(summary_h2$file,"_"),"[[",3)
summary_h2$stat<-"h2"
# add flag for TBV adjustment
summary_h2$adj<-""
summary_h2$adj[grep("TBV",summary_h2$file)]<-"adjTBV"

# subset of significan AIs - pairs
h2AI_sig<-subset(summary_h2,(pval<0.05))

# Combine both:
summary<-summary_h2
summary$sig<-"ns"
summary$sig[summary$pval<0.05]<-"p<0.05"
summary$sig<-factor(summary$sig, c("ns","p<0.05","p<0.01","p<0.0002"))#levels=c("ns","p<0.05")) # 

#----------------------------------------------------------------------

#----------------------------------------------------------------------
# plots
#----------------------------------------------------------------------
mytheme<-theme_bw() + theme(panel.spacing = unit(0, "lines"), 
                            strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                            axis.title.x=element_blank(), axis.ticks.x =element_blank(), # axis.text.x =element_blank(),
                            axis.text.x = element_text(angle = 90, size=16,colour="black",hjust=1,vjust=.5),
                            axis.title.y=element_text(size=18),axis.text.y=element_text(size=16,colour="black"),legend.text=element_text(size=16)) 
mytheme2<-theme_bw() + theme(panel.spacing = unit(0, "lines"), 
                            strip.background = element_rect(fill="white"), strip.text = element_text(size=16), 
                            # axis.title.x=element_blank(), axis.ticks.x =element_blank(), # axis.text.x =element_blank(),
                            # axis.text.x = element_text(size=16,colour="black",hjust=1,vjust=.5),
                            # axis.title.y=element_text(size=16),axis.text.y=element_text(size=16,colour="black"),
                            title=element_text(size=16),
                            axis.title=element_text(size=16),axis.text=element_text(size=14,colour="black"),
                            legend.text=element_text(size=16), legend.title =element_text(size=16) ) 

cols=c("#d95f02","#7570b3","#1f78b4","#33a02c","#000000") # '#1b9e77','#d95f02','#7570b3 # "#a6cee3"  ,"#7570b3"
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#----------------------------------------------------------------------
# meas_labels<-unlist(c(bquote(''*h^2*'(AI native)'),
#                       bquote(''*h^2*'(AI symmetric)'),
#                       bquote(''*h^2*'(L)'),bquote(''*h^2*'(R)'),bquote(''*rho*'(L,R)')))
# 
# meas_labels2<-c('AI native','AI symmetric','L','R')

# need to fix: header, colors ... spacing between regions
h2_plot<- ggplot(data=summary,aes(x=Dist,y=h2,color=phenotype)) +
  geom_pointrange(aes(alpha=sig,shape=sig,ymin=h2-se, ymax=h2+se),position = position_dodge(width=0.5),size=1,stat="identity") + 
  facet_grid(adj~Asy,scales="free_x",space = "free_x") + #facet_grid(.~Asy,scales="free_y") +
  mytheme + 
  ylab(bquote('Estimate ('*h^2*')')) +
  scale_color_manual(values=cols,name = "")  +
  scale_alpha_discrete(range = c(0.4, 1), name="") + 
  scale_shape_discrete(name="") + coord_cartesian(ylim=c(0,1))  + theme(legend.position="bottom")
# save output
ggsave(h2_plot,file=paste(out_dir,"Skews.png",sep=""),width=10)
write.csv(summary,file=paste(out_dir,"Skews.csv",sep=""),quote=FALSE,row.names=FALSE)
#----------------------------------------------------------------------

