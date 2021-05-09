# create input data for heritability analyses using FS derived asymmetry phenotypes by Xiangzhen Kong
library(pander)
library(ggplot2)
library(car)
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(knitr)
library("GGally")

opts_chunk$set(include=TRUE, warnings=FALSE, echo=FALSE, results = "asis", tidy=TRUE, width=50, fig.width=8, fig.height=6)

panderOptions('knitr.auto.asis',FALSE)

# Define directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {
  dir="P://workspaces/"
} else {
    dir="/data/workspaces/lag/workspaces/"
    }

primary_dir=paste(dir,"lg-ukbiobank/analysis/xiangzhen/brainTorque/torque_ukb_N40681/out_2/",sep="") # where pheno data is stored, also update mfiles far below

working_dir=paste(dir,"lg-ukbiobank/analysis/xiangzhen/brainTorque/torque_ukb_N40681/GWAS_Dick/GWAS_out2/",sep="")
dir.create(file.path(working_dir), showWarnings = FALSE)

out_plots_dir=paste(working_dir,"plots/")
dir.create(file.path(out_plots_dir), showWarnings = FALSE)

setwd(working_dir)

# get functions to calculate AIs, some histogram wrap, etc
source(paste(dir,"lg-ukbiobank/analysis/amaia/phenotypes/AIfunctions.R",sep=""))

###  -----1-----   ###

# read data containing all the demographic info and not-transformed phenotypes
# created when running: ukb9246_ukb10785_selectCols_subsetSamples.R
#amaia_working_dir=paste(dir,"lg-ukbiobank/working_data/amaia/pheno_files/FSphenos_xiakon/",sep="")
#data<-read.csv(
#  paste(working_dir,"../samples/demographic/",
#        "ukb40473_sqc_fam_selectedPhenotypes_imagingSamples.csv",
#        sep=""))
data<-read.csv(
  paste(working_dir,"../demographic/",
        "sqc_sample_ukb40473.csv",
        sep=""))

# read sample QC info, and add column
#qc<-read.table(paste(working_dir,"../samples/ukb9246_imagingT1_N40681_samples.list", sep=""))
#qc$QC<-1
#data_qc<-merge(data,qc[,-2],by.x="f.eid",by.y="V1",all.x=TRUE)
data_qc<-data
data_qc$QC<-1
table(data_qc$QC) # 8400 that pass genetic QC, and relatedness
data<-data_qc; rm(data_qc)

###  -----2-----   ###

#------------------------------------------------------------------------------------
#' select covariates and useful info cols for filtering

#' Age at assesment, instance 2 --> age at imaging visit
age1<-colnames(data)[grep("^Age_when",colnames(data))] # this is average age across three visits
age<-"f.21003.2.0" # age at imaging visit
print("is.na age")
table(is.na(data[,age1])) # no missing data
table(is.na(data[,age])) # missing data for several samples... 
# if we include this variable in the lm, will remove ~ 750 samples

# calculate age from birth date and date attending assesment center (imaging visit): f.53.2.0
# birth year:  f.34.0.0
# birth month: f.52.0.0
data$birth_y<-as.character(data$f.34.0.0)
data$birth_m<-as.character(data$f.52.0.0)
data$birth_m[nchar(data$f.52.0.0)==1]<-paste("0",data$f.52.0.0[nchar(data$f.52.0.0)==1],sep="")
data$birth_ym<-as.Date(paste(data$birth_y,data$birth_m,"01",sep="-"),format="%Y-%m-%d")
data$date_imaging_ym<-as.Date(data$f.53.2.0,format="%Y-%m-%d")
data$date_imaging_y<-as.numeric(sapply(strsplit(as.character(data$f.53.2.0),"-"),"[[",1))
# calculate at age imaging, taking into account year and month
data$age_imaging<-round(as.numeric(data$date_imaging_ym-data$birth_ym)/365.25,digits=0)
data$age_imaging1<-round(as.numeric(data$date_imaging_ym-data$birth_ym)/365.25,digits=1)
data$age_imaging_f<-floor(as.numeric(data$date_imaging_ym-data$birth_ym)/365.25)
data$age_imaging_y<-data$date_imaging_y-as.numeric(data$birth_y)
plot(data$age_imaging,data[,age])
plot(data$age_imaging_y,data[,age])
plot(data$age_imaging_f,data[,age])
plot(data$age_imaging_y,data$age_imaging)
plot(data$age_imaging_y,data$age_imaging_f)
plot(data$age_imaging_y,data$age_imaging1)
# compare to reported age
head(data[which(data$age_imaging-data[,age]!=0),c(age,"f.53.2.0","birth_y","birth_m","date_imaging_ym","date_imaging_y","age_imaging","age_imaging_y","age_imaging_f","age_imaging1")])

# go with age_imaging1, more precision...
data$age<-data[,"age_imaging1"]
#age<-"f.21003.2.0"  # ------------------------- use "f.21003.2.0" for imaging visiting age; use "age" for estimated one
age<-"age"  # ------------------------- use "f.21003.2.0" for imaging visiting age; use "age" for estimated one

#' Assesment center
#assesmentC<-colnames(data)[grep("UK_Biobank_assessment",colnames(data))] # missing, need to add!
assesmentC<-colnames(data)[grep("f.54.2.0",colnames(data))] # Corrected!

# T1 quality
#t1cnr<-colnames(data)[grep("Inverted_contrast.to.noise_ratio_in_T1",colnames(data))]
#t1snr<-colnames(data)[grep("Inverted_signal.to.noise_ratio_in_T1",colnames(data))]
t1cnr<-'f.25735.2.0'
t1snr<-'f.25734.2.0'

#' Sex  
sex<-"f.31.0.0"
boxplot(data[,age]~data[,sex]) # males are older than females on average
print("is.na sex")
pander(table(is.na(data[,sex]))) 
#' Some missing data, same as people without Genetic info.  
#' Select Genotyping array and batch.  
#array<-colnames(data)[grep("genotyping.array",colnames(data))]
batch<-colnames(data)[grep("f.22000.0.0",colnames(data))] # genotyping batch
data$array<-(data[,batch]>0)*1
array<-'array'

# scanner position parameters
# scanner_position<-colnames(data)[grep("brain_position",colnames(data))]
scanner_position<-colnames(data)[grep("25756.2.0|25757.2.0|25758.2.0",colnames(data))]

# total brain volume
#tBV<-"Volume_of_brain._grey.white_matter"#colnames(data)[grep("Volume_of_brain._grey.white",colnames(data))]
#tBV<-colnames(data)[grep("Volume_of_brain._grey.white",colnames(data))]
tBV<-'f.25010.2.0'
# motion covs
#motion<-colnames(data)[grep("rfMRI_head_motion",colnames(data))]
motion<-'f.25741.2.0'

#' ## Define (not imaging) phenotypes:  
# handedness<-colnames(data)[grep("^handedness_",tolower(colnames(data)))] # this was mean across instances were handedness was assessed
handedness<-colnames(data)[grep("^f.1707.",tolower(colnames(data)))]
handedness<-handedness[1]
# define congruent handedness, right or left only
table(data[,handedness]) #coding is as factor: 1=rigth-handed; 2=left-handed
data$handedness<-NA
data$handedness[data[,handedness]==1]<-0
data$handedness[data[,handedness]==2]<-1 # left-handers as cases
table(data$handedness)
#' Left handers (consistent): `sum(data$handedness==1,na.rm=TRUE)` out of `NROW(data)`

# compute Zage2, as: ((age-meanAge)/stdAge)2
meanAge<-mean(data[,age],na.rm = TRUE) # ] 62.73657
stdAge<-sqrt(var(data[,age],na.rm = TRUE)) # ] 62.73657
data$zage2<-((data[,age]-meanAge)/stdAge)^2 # age2: (age-meanAge)^2
rm(meanAge)
# check
plot(data$age,data$zage2)

place_birth1 = 'f.129.0.0'
place_birth2 = 'f.130.0.0'
place_home1 = 'f.20074.0.0'
place_home2 = 'f.20075.0.0'
weight = 'f.21002.2.0'
height = 'f.12144.2.0'

#' covariates to keep
N_PCs = 10#10 40  # Edit! --- Number of PCs to control for. ---
covs2keep<-c(assesmentC,age,"zage2",sex,paste("f.22009.0.",1:N_PCs,sep=""),array,batch,scanner_position, tBV, motion, t1cnr, t1snr)#, place_birth1, place_birth2, place_home1, place_home2, height, weight)

#' Define other phenotypes of interest
phenos2keep<-c("handedness") # add unique to avoid having duplicated columns!
#------------------------------------------------------------------------------------

# subset columns of interest
data2<-data[,c("f.eid","f.eid","QC",covs2keep,phenos2keep)] # ; rm(data)
colnames(data2)<-c("ID1","ID2","QC",
                   "assessment_centre","age","zage2","sex",
                   paste("PC",1:N_PCs,sep=""),"array","batch",
                   "Scanner lateral (X) brain position",
                   "Scanner transverse (Y) brain position",
                   "Scanner longitudinal (Z) brain position",
                   "totalBV",#"totalBVadjHeadSize",
                   "motion", 't1cnr', 't1snr',
                   "handedness"#,
                   #'place_birth1', 'place_birth2', 'place_home1', 'place_home2', 'height', 'weight'
                   )
# fix names
colnames(data2)<-gsub(" ","_",gsub("\\(|\\)|\\.","",colnames(data2)))

# redefine covariate variables, just in case
#age<-"age"
scanner_position<-colnames(data2)[grep("Scanner_",colnames(data2))]
motion<-colnames(data2)[grep("motion",colnames(data2))]
#tBV<-colnames(data2)[grep("grey*.*white",colnames(data2))]
bin_covs2<-c("sex","assessment_centre","array","batch")
quant_covs2<-c("age","zage2",paste("PC",1:N_PCs,sep=""),scanner_position,"totalBV", 't1cnr', 't1snr')#, 
               #'place_birth1', 'place_birth2', 'place_home1', 'place_home2', 'height', 'weight')

###  -----3-----   ###

#------------------------------------------------------------------------------------
#' Visually inspect covariates to see whether there are outliers # if so, should blank?!
quant_covs_hist_out<-lapply(quant_covs2,function(x){histo_out_AIs(data2=data2,x,thr=8)})
# save plot
hists<-sapply(quant_covs_hist_out,"[[",2)
pdf(paste(out_plots_dir,"quant_covariates_outliers.pdf",sep=""),onefile = TRUE)
for (h in 1:length(hists)){ grid.arrange(hists[[h]]) }
dev.off()
rm(hists)
#' there are two very clear outliers in scanner_position_x; blank them # more??? **
w<-which(data2[,scanner_position[1]]<(-20))
if (length(w)>0) {data2[w,scanner_position[1]]<-NA}
rm(w)

###  -----4-----   ###

#------------------------------------------------------------------------------------
#' Define order of covariates to include within lm
#' Drop batch, too many levels # ,"batch" ; check assessment centre
#' drop assessment center as well# "assessment_centre",

# # witout adjusting for TBV

covs_order<-c("age","zage2","sex",paste("PC",1:N_PCs,sep=""),"array",scanner_position, 'assessment_centre', 't1cnr', 't1snr')#,'place_birth1', 'place_birth2', 'place_home1', 'place_home2', 'height', 'weight')  #Edit!!!
covs_name="_noBioCovs_ZAge"#<---GWAS_out


# defines covs, # adjusting for totalBV as well
#covs_order<-c(covs_order,"totalBV")
#covs_name="_noBioCovs_ZAge_TBV"
#covs_name="_noBioCovs_noAssessmentC_TBV"
#covs_name="_noBioCovs_ZAgeTBV" # Edit!!!
#covs_name="_noBioCovs_ZAgePC40_TBV" # <--- GWAS_out3

# reorder covariates, only if present in covs_order, i.e. included in LM
quant_covs2<-quant_covs2[quant_covs2 %in% covs_order]
bin_covs2<-bin_covs2[bin_covs2 %in% covs_order]

###  -----5-----   ###

#------------------------------------------------------------------------------------
#' phenotypes derived by Xiangzhen

# identify files:
mfiles<-paste(primary_dir,list.files(primary_dir,pattern="SkewsXY*"),sep="")  #Edit!!!

mfiles<-mfiles[grep("NumVert",mfiles,invert=TRUE)]

for (f in mfiles) {
  # extract info from file name
  f_info<-sapply(strsplit(gsub(primary_dir,"",f),"_"),"[[",1)
  # read data
  dat<-read.csv(f,header=TRUE)
  w<-which(colnames(dat)=="SID")
  phenos<-colnames(dat)[-w]
  
  # create directories for output data
  dir.create(file.path(working_dir,f_info), showWarnings = FALSE)
  dir.create(file.path(paste(working_dir,f_info,sep="/")), showWarnings = FALSE)
  dir.create(file.path(paste(working_dir,f_info,sep="/"),"lm"), showWarnings = FALSE)
  dir.create(file.path(paste(working_dir,f_info,sep="/"),"lm_plots"), showWarnings = FALSE)
  output_dir=paste(working_dir,f_info,sep="/")
  
  # combine data with covariates
  d_asy<-merge(data2,dat,by.x="ID1",by.y="SID",all.y=TRUE); rm(dat)
  table(d_asy$QC) # N = 7301
  data_s<-subset(d_asy,QC==1); rm(d_asy)
  
  #' Convert binary covariates into factors
  data_s[,bin_covs2]<-apply(data_s[,bin_covs2],2,function(x) as.factor(x))
  data_s<-subset(data_s,!is.na(data_s$sex))
  table(data_s$sex)
  # pander(str(data_s[,phenos]))
  data_s$handedness<-factor(data_s$handedness,levels=c(0,1))
  table(data_s$handedness)
  #' For every phenotype in data_s, blank outliers: if +-4SD -------------------------------------------outliers
  t=4
  out_blank<-lapply(data_s[,phenos], function(x,thr=t){
    if (class(x)!="factor"){
      m<-mean(x,na.rm=TRUE)
      sd<-sd(x,na.rm=TRUE)
      out<-which(x>(m+thr*sd)|x<(m-thr*sd))
      if(length(out)>0) {x[out]<-NA}
    }
    return(x)
  })
  out_blank<-do.call("cbind",out_blank)
  counts_blanked<-apply(out_blank,2,function(x) sum(is.na(x)))
  pander(counts_blanked)
  #' Maximum number of outliers blanked using a threshold of +- `t` SD from the mean: `max(counts_blanked)`
  write.csv(counts_blanked,file=paste(output_dir,"/",f_info,"_counts_blanked_4SD",covs_name,".csv",sep=""),row.names = TRUE,quote=FALSE)
  rm(out_blank,counts_blanked)
  
  for (p in phenos) {
    print(paste("Dependent variable: ",p,"\n  ",sep=""))
    # for each phenotype, first check if output exists and only proceed if it does not
    lm_file=paste(output_dir,"/lm/",p,"_lm",covs_name,".txt",sep="")
    # also, check that it's not all 0's or NAs
    if (sum(!(is.na(data_s[,p])|data_s[,p]=="0")) >0) {
    # if (file.exists(lm_file)==FALSE) {
      # define name of column for residuals
      p2<-paste("residuals_",p,sep="")
      data_s[,p2]<-data_s[,p] # and save data there
      #----------------------
      print('Detect and blank outliers"\n  ')
      m<-mean(data_s[,p],na.rm=TRUE)
      sd<-sd(data_s[,p],na.rm=TRUE)
      out<-which( data_s[,p] > m+t*sd |data_s[,p] < m-t*sd )
      print('Remove outliers: ')
      print(length(out))
      print('')
      data_s[out,p2]<-NA
      rm(m,sd,out)
      #----------------------
      #print('Regress out covariates and save residuals"\n  ')
      print(paste('Run a linear model to get residuals for ',p,'\n  ',sep=""))
      model_formula<-paste(p2,"~",paste(c(covs_order),collapse=" + "))
      if (class(data_s[,p2])=="numeric" | class(data_s[,p])=="integer") {
        lm<-do.call("lm",list (as.formula(model_formula),data=data_s))
      } else if (class(data_s[,p2])=="factor") {
        lm<-do.call("glm",list (as.formula(model_formula),data=data_s,family=binomial(link='logit')))
      }
      pander(anova(lm),round=1000)
      #----------
      print('Fit check\n  ')
      # par(mfrow=c(2,2))
      pdf(file=paste(output_dir,"/lm_plots/",p,"_lm_plots",covs_name,".pdf",sep=""))
      plot(lm)
      dev.off()
      par(mfrow=c(1,1))
      #----------
      # save model + summary model estimates in txt file as well
      df<-as.data.frame(anova(lm)); df$variable<-rownames(df)
      df<-df[,c(NCOL(df),1:(NCOL(df)-1))]
      df2<-as.data.frame(summary(lm)$coefficients); df2$variable<-rownames(df2)
      df2<-df2[,c(NCOL(df2),1:(NCOL(df2)-1))]
      write.table(df,file=paste(output_dir,"/lm/",p,"_lm",covs_name,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
      write.table(df2,file=paste(output_dir,"/lm/",p,"_coefficients_lm",covs_name,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)
      rm(df,df2)
      #----------
      tmp_df<-data_s[,c("ID1","ID2",p2,bin_covs2,quant_covs2)]
      covNA<-names(which(apply(tmp_df,1,function(x) sum(is.na(x)))>0))
      w<-match(covNA,rownames(tmp_df)); rm(covNA)
      length(w)
      if (length(w)>0) { tmp_df<-tmp_df[-w,] }
      tmp_df$predicted=predict(lm)
      tmp_df$residuals=residuals(lm)
      #----------
      m<-mean(tmp_df$residuals,na.rm = TRUE)
      sd<-sd(tmp_df$residuals,na.rm = TRUE)
      res_out<-which(tmp_df$residuals > m+t*sd | tmp_df$residuals < m-t*sd )
      print('Number of outliers for residuals, deviating more than 4SD from mean:')
      print(length(res_out)) # do not blank them again!
      # if (length(res_out) > 0) {
      #   tmp_df$residuals[res_out]<-NA
      # }
      rm(m,sd,res_out)
      # combine residuals with data_s
      data_s[,p2]<-NA
      if (length(w)>0) { data_s[-w,p2]<-tmp_df$residuals; data_s[w,p2]<-NA
      } else { data_s[,p2]<-tmp_df$residuals }
      rm(w)
      #----------
      # some more plots
      if (class(data_s[,p])!="factor") {
        print('More model plots\n  ')
        # plot phenotype in relation to binary covariates, boxplots
        print(paste('Plot ',p,' in relation to binary covariates, boxplots\n  ',sep=""))
        for (cov in bin_covs2) {
          if (class(data_s[,p])=="numeric"| class(data_s[,p])=="integer") {
            plot<-ggplot(data=data_s,aes_string(x=cov,y=p)) + geom_violin() +geom_boxplot(width=0.1)+ theme_bw()
          } else if (class(data_s[,p])=="factor"){
            plot<-ggplot(data=subset(data_s[-which(is.na(data_s[,p])),],!is.na(sex))) +
              geom_bar(aes_string(x=cov,fill=p),position="fill") +
              scale_fill_grey(start=.3,end=.7) + theme_bw()
          }
          assign(paste(p,cov,"plot",sep="_"),plot) ; rm(plot)
        }
        # pander(anova(lm),round=100)
        ncol=round(length(bin_covs2)/2)
        list_p<-lapply(ls(pattern=paste(p,"*.*plot",sep="")),get)
        g<-grid.arrange(arrangeGrob(grobs=list_p))
        ggsave(file=paste(output_dir,"/lm_plots/",p,"_plot_binary_covs",covs_name,".pdf",sep=""),g,width=10,height=10)
        rm(cov,ncol,list_p,g)
        rm(list=ls(pattern="plot"))
        # plot phenotype in relation to quantitative covariates, scatterplot
        print(paste('Plot ',p,' in relation to quantitative covariates\nScatterplot with lm and loess smooth',sep=""))
        for (cov in quant_covs2) {
          if (class(data_s[,p])=="numeric"| class(data_s[,p])=="integer") {
            plot<-ggplot(data=data_s, aes_string(x=cov,y=p)) + geom_point() +
              geom_smooth(span=0.3,method='loess') + geom_smooth(span=0.3,method='lm',colour='red') +
              theme_bw()
          } else if (class(data_s[,p])=="factor"){
            plot<-ggplot(data=data_s, aes_string(x=p,y=cov)) +
              geom_boxplot(varwidth = TRUE) +
              theme_bw()
          }
          assign(paste(p,cov,"plot",sep="_"),plot) ; rm(plot)
        }
        ncol=round(length(quant_covs2)/2)
        list_p<-lapply(paste(p,quant_covs2,"plot",sep="_"),get)
        g<-grid.arrange(arrangeGrob(grobs=list_p))
        ggsave(file=paste(output_dir,"/lm_plots/",p,"_plot_quant_covs",covs_name,".pdf",sep=""),g,width=10,height=10)
        rm(cov,ncol,list_p,g)
        rm(list=ls(pattern="plot"))
        rm(list=ls(pattern="^p_"))
        if (class(data_s[,p])=="numeric"| class(data_s[,p])=="integer") {
          # plot histogram of original phenotype and residuals
          p_hist<-ggplot(data=data_s, aes_string(x=p)) + geom_histogram() + geom_vline(xintercept = mean(data_s[,p])-sd(tmp_df[,p2])*t) + geom_vline(xintercept = mean(data_s[,p])+sd(data_s[,p])*t) + theme_bw()
          p_rhist<-ggplot(data=data_s, aes_string(x=p2)) + geom_histogram() +  theme_bw()
        } else if (class(data_s[,p])=="factor") {
          p_hist<-ggplot(data=data_s, aes_string(x=p)) + geom_bar() + theme_bw()
          p_rhist<-ggplot(data=data_s, aes_string(x=p2)) + geom_histogram() + theme_bw()
        }
        list_p<-lapply(ls(pattern="^p_"),get)
        g<-grid.arrange(arrangeGrob(grobs=list_p))
        ggsave(file=paste(output_dir,"/lm_plots/",p,"_histograms",covs_name,".pdf",sep=""),g,width=10,height=10)
        rm(list=ls(pattern="^p_"))
        #
        # print('Get residuals and save for downstream analyses\n')
        # file=paste(working_dir,"/pheno_files/genetic_v2/ukb9246_ukb10785_imaging",covs_name,"_residuals.",p,"_wHeader.table",sep="")
        # write.table(tmp_df[,c("ID1","ID2","residuals")],file=file,sep="\t",quote=FALSE,row.names = FALSE,col.names=TRUE)
        # # clean
        rm(lm,tmp_df,file,lm_file)
      }

      }
  }
 
  
  # save list of phenotypes
  name<-gsub(".out","_residuals4genetics",gsub(primary_dir,"",f))
  
  # save phenotypes
  write.table(data_s,paste(working_dir,"/",name,covs_name,".table",sep=""),col.names = TRUE,row.names = FALSE,quote=FALSE)
  # save list of residuals
  residuals_list<-colnames(data_s)[grep("residuals",colnames(data_s))]
  write.table(residuals_list,paste(working_dir,"/",name,covs_name,".list",sep=""),col.names = FALSE,row.names = FALSE,quote=FALSE)
  # clean
  rm(name,residuals_list)
  rm(data_s)

}

###  -----6-----   ###

rm(list=ls())
gc()

