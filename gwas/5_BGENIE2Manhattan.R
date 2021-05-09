if("qqman" %in% rownames(installed.packages()) == FALSE) {install.packages("qqman")}
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {install.packages("gridExtra")}
library(qqman)
library(gridExtra)
library("data.table")
options(stringsAsFactors = FALSE) #,digits=22

# set number of digits to 22, default is 7, but this affects calculation of lambda, which always becomes 1
#---------------------------------------------------#
# args=commandArgs()
# file=args[1]
# p_col=args[2]
#---------------------------------------------------#
# It somehow fails (due to gzfile I think) if I try to run it directly from the terminal..., so run a loop on all output
setwd("/home/xiakon/MPI_workspace/lg-ukbiobank/analysis/xiangzhen/brainTorque/torque_ukb_N40681/GWAS_Dick/GWAS_out_withplaceheightweight/output") # 
# create ./plots ./clean directories if they do not exist
if (!file.exists("plots")){
  dir.create(file.path("plots"))
}

if (!file.exists("clean")){
  dir.create(file.path("clean"))
}

# read filtered snps based on HWE<1e-07 and INFO<0.7 (INFO from QCtool - total sample, not BGENIE!)
qc_dir="/data/workspaces/lag/workspaces/lg-ukbiobank/working_data/genetic_data/subset_imagingT1_40k/v2_white_ancestry/with_rel_filter/VQC/"

# Parse arguments:
## chr, phenotype
#------------------------------------
# Get command line arguments
args = commandArgs(trailingOnly=TRUE)

phenos2test<-c('SkewsXY', 'SkewsXZ')


for (phenotype in phenos2test){

  args<-c("all",phenotype,
        "GWAS_Dick",
        "imagingT1_N40681")

  chr=args[1] #either all, or a specific chr name
  chr="all"
  pheno=gsub("residuals_","",args[2]) # phenotype name
  pheno_root=args[3] # ukb phenotype batch names
  subset_name=args[4] # subset of data that will be included,, in case not samples are in; should reflect sample_file info

  #------------------------------------
  # define thresholds
  maf_thr=0.01#0.01#0.001 # args[6]
  info_thr=0.7 # args[7]
  hwe_thr=1e-7
  #------------------------------------
  rm(args)


  key=""    #"TBV"  # set to "", if not adjusted for TBV
  # define phenotype name
  name=gsub("_"," ",gsub("residuals_","",pheno))
  if (key!=""){
    name=paste(name,key,sep="")
  }
  #------------------------------------

  print("Check if Manhattan plot (or some other output file) exists, only proceed if it does not")

  if (file.exists(paste("plots/Manhattan_",gsub(" ","_",name),"_CHR",chr,".png",sep=""))==FALSE){

  if (chr!="all") {c=paste(chr,"\\.",sep="")} else {c=".*."} # just one chr!
  # define files to parse, will be either one (if one chr, or all)
  pattern1=paste(pheno_root,subset_name,key,"bgenie","chr",sep="_")
  pattern4files<-paste(pattern1,c,sep="")
  files<-list.files(pattern=pattern4files)
  rm(pattern4files,c)

  cat('Frequency filter threshold:',maf_thr,'\n')
  cat('Info filter threshold: ',info_thr,'\n')
  cat('HWE p-value filter threshold:',hwe_thr,'\n')

  # function to parse files, read one at a time, and subset only cols for phenotype of interest
  for (file in files) {
    # define chromosome
	cat(file,'\n')
    chrom=gsub(".out.gz","",gsub(pattern1,"",file))
    #-----------------------------------
    # get snp QC file, to filter out: HWE<1e-7 or INFO<0.7
    qc_f<-paste(qc_dir,list.files(qc_dir,pattern=paste("_chr",chrom,".snpstats_mfi_hrc.txt",sep="")),sep="")
    #if (chrom!="X"){
    #  qc_f<-paste(qc_dir,list.files(qc_dir,pattern=paste("_chr",chrom,"_v3*.*_snp-stats_mfi_hrc.txt",sep="")),sep="")
    #} else {
    #  qc_f<-paste(qc_dir,list.files(qc_dir,pattern=paste("_chr",chrom,"_v3*.*_snp-stats_infer_ploidy_mfi_hrc.txt",sep="")),sep="")
    #}
    snp_qc<-fread(qc_f,header=TRUE)
    if (chrom=="X"){
      snp_qc$HW_exact_p_value.subset.QCtool<-snp_qc$HW_females_exact_pvalue.subset.QCtool
      snp_qc$HW_lrt_p_value.subset.QCtool<-snp_qc$HW_females_lrt_pvalue.subset.QCtool
    }
  # qc_f<-paste(qc_dir,list.files(qc_dir,pattern=paste("_chr",chrom,"_v3*.*_snpsFiltered_hwe1eMin7_info0.7.txt",sep="")),sep="")
  # if (chrom!="X"){
  #   # colnames(snp_qc)<-c("chr_pos","rsid","hwe","info.qctool","missingness","total")
  #   h<-read.table(paste(qc_dir,list.files(qc_dir,pattern=paste("chrA","_*.*header",sep="")),sep=""))
  #   colnames(snp_qc)<-as.character(h[1,])
  # } else {
  #   h<-read.table(paste(qc_dir,list.files(qc_dir,pattern=paste("chr",chrom,"_*.*header",sep="")),sep=""))
  #   colnames(snp_qc)<-as.character(h[1,])
  # }
  # rm(qc_f,h)
  #-----------------------------------
  # start processing file
  # cat(paste("Opening conection for chromosome ",chrom,"\nFile: ",file,'\n',sep=""))
    zz=gzfile(file,'rt')  
  # get header first, and define columns of interest
    cat('Reading header of file...\n')
    header=read.table(zz,header=T,skipNul=TRUE,nrows=1)
  
  # define general cols
    snp_cols<-grep("residuals",colnames(header),invert=TRUE)
    pheno_cols<-grep(gsub("_totalBV","",pheno),colnames(header))
  
  # read file
    cat(paste('Reading file...\n'))
    cat(as.character(Sys.time()),'\n')
    d0=fread(paste('zcat ',file,sep=""),header=T,select=c(snp_cols,pheno_cols))
    cat(as.character(Sys.time()),'\n')
    cat('Number of SNPs in choromsome ',chrom,':',NROW(d0),'\n')
    #  combine with qctool snp qc info
    print(colnames(d0))
    print(colnames(snp_qc))
    d0=merge(d0,snp_qc,
             by.x=c("rsid","chr","pos","a_0","a_1"),
             by.y=c("RS_ID.UKB","CHR","Position","A1.UKB","A2.UKB"),
             all.x=TRUE)
    table(d0$filter) # number of SNPs that should be filtered out...
  # check filters:
  # maf
    cat('Number of variants with af<',maf_thr,'\n')
    cat(table(d0$af<=maf_thr))
    cat('Number of variants with af>',(1-maf_thr),'\n')
    cat(table(d0$af>=(1-maf_thr)))
    cat('Number of filtered out based on frequency (TRUE) for chromosome',chrom,'\n')
    cat(table(d0$af<=maf_thr|d0$af>=(1-maf_thr)))
  # info
    cat('Number of variants with info <',info_thr,'\n')
    cat(table(d0$info<=info_thr),'\n')
    cat(table(d0$info.imaging.QCtool<=info_thr),'\n')
    cat(table(d0$INFO.UKB<=info_thr),'\n')
    #table(UKB.MFI=d0$INFO.UKB<=info_thr,subset=d0$info.imaging.QCtool<=info_thr)
  # subset on maf
    d1<-subset(d0,(af>=maf_thr&af<=(1-maf_thr)))
    cat('Number of SNPs in choromsome ',chrom,' after filtering:',NROW(d1),'\n')
    cat('Number of variants with HWE pval <',hwe_thr,'\n')
    cat(table(d1$hwe<hwe_thr))
    cat('Number of variants with INFO.qctool  <',info_thr,'\n')
    cat(table(d1$info.qctool<info_thr))
  # double check INFO field from qctool, and hwe
    cat('Filter out variants with UKB.INFO<', info_thr,' or HWE pval <', hwe_thr)
    d2<-subset(d1,INFO.UKB>info_thr&HW_exact_p_value.subset.QCtool>hwe_thr) # only 9367 markers in chrX
    cat('Number of variants filtered out for chromosome ', chrom,':')
    NROW(d1)-NROW(d2) # variants filtered out because of hwe / INFO
  # some extra checks for chr X
  # some extra checks for chr X
    if (chrom=="X"){
      d2<-subset(d1,INFO.UKB>info_thr&HW_exact_p_value.subset.QCtool>hwe_thr) 
      d2<-d2[,colnames(d),with=FALSE]
      w=which(colnames(d2)=="female")
      if(length(w)>0) {
      d2<-d2[,-w,with=FALSE]
      }
      rm(w)
    }
    cat('Number of variants filtered out for chromosome ', chrom,':')
  # save 1/d2 into all chromosomes
    if (!exists("d")) {d<-d2} else {d<-rbind(d,d2)}
  
  # clean
    rm(chrom,zz,header,snp_cols,pheno_cols)
    rm(d0,d1,d2)
  } ; rm(file,files)

# select columns to keep, otherwise there are too many
  cols<-colnames(d)[grep(".QCtool|.UKB|.BGENIE|.HRC",colnames(d),invert=TRUE)]
  cols2<-colnames(d)[grep(".QCtool|.UKB",colnames(d),invert=FALSE)]
  cols2<-cols2[grep("hw_exact_p_value.|^info\\.|maf\\.",tolower(cols2))]
#
  data<-d[,c(cols,cols2),with=FALSE]
  rm(cols,cols2)
  rm(snp_qc,d)
# check SNPs per chr
  dim(data)
  table(data$chr)
#-----------------------------------------------------------

# rename columns
  colnames(data)<-gsub("^_","",gsub(pheno,"",gsub("residuals_","",colnames(data))))
# get p from -logP or beta!
  data$P<-10^(-data$`-log10p`)
  data$P_z<-2*pnorm(-abs(data$beta/data$se))
  plot(data$P,data$P_z)
  p_col="P"
# subset
  head(data)
  data1<-subset(data,af>0.001&af<0.999&info>info_thr)
# data1<-subset(data,A1FREQ>0.01&A1FREQ<0.99&INFO>0.7)
  data_excluded<-subset(data,!(af>0.001&af<0.999&info>0.7))
# subset to plot in Manhattan
  data2<-subset(data1,P<0.01)
  data_sig<-subset(data2,P<5e-07)
  write.csv(data_sig,file=paste("clean/",gsub(" |\n","_",name),"_CHR",chr,"_Pmin7.csv",sep=""),row.names = FALSE)
  write.table(data1,file=paste("clean/",gsub(" |\n","_",name),"_CHR",chr,"_",hwe_thr,"HWEp_",info_thr,"INFO_",maf_thr,"MAF.txt",sep=""),row.names = FALSE,col.names=TRUE,quote=FALSE)

#---------------------------------------------------#
# Plots
#---------------------------------------------------#
# For p-values, calculate chi-squared statistic

  chisq <- qchisq(1-data$P,1)
  mchisq<-summary(chisq)[c("Mean","Median","Max.")]
  chisq1 <- qchisq(1-data1$P,1)
  mchisq1<-summary(chisq1)[c("Mean","Median","Max.")]

# Calculate lambda gc (??gc)
  print("Calculate lambda")
  lambda <- round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),digits=5)
  lambda1 <- round(median(chisq1,na.rm=TRUE)/qchisq(0.5,1),digits=5)
  print(lambda)
  print(lambda1)
# save lambda values into file
  t<-cbind(name=file,
         lambda_all=lambda,
         mean_chisq_all=mchisq[1],median_chisq_all=mchisq[2],max_chisq_all=mchisq[3],
         lambda_clean=lambda1,
         mean_chisq_clean=mchisq1[1],median_chisq_clean=mchisq1[2],max_chisq_clean=mchisq1[3]
         )
  write.csv(t,file=paste("clean/",gsub(" |\n","_",name),"_CHR",chr,"_lambdas_chisqStats.csv",sep=""),row.names = FALSE)
  rm(t)

  pmin=ceiling(abs(log10(min(data2$P))))
  data2$chrom<-data2$chr
  data2$chrom[data2$chr=="X"]<-23
  data2$chrom<-as.numeric(data2$chrom)
  print("Plot\n")
  png(paste("plots/Manhattan_",gsub(" |\n","_",name),"_CHR",chr,".png",sep=""),width=1000)#, height=1000,
  # pdf(paste("plots/Manhattan_",gsub(" |\n","_",name),"_CHR",chr,".pdf",sep=""),width=15) #, height=1000,
  manhattan(data2,chr="chrom",bp="pos",p="P",snp="rsid",ylim=c(2,pmin),main=name,
          suggestiveline= (-log10(5e-8)),genomewideline= (-log10(1e-11)))

  dev.off()

  # pdf(paste("plots/QQplot_",gsub(" |\n","_",name),"_CHR",chr,".pdf",sep=""))

  png(paste("plots/QQplot_",gsub(" |\n","_",name),"_CHR",chr,".png",sep=""))
  qq(data1$P,main=paste("Q-Q plot of GWAS p-values\n",name,sep=""))
  text(x=3,y=0.5,bquote(lambda==  .(lambda1)))
  dev.off()
  print("Finished, hopefully without errors")
  # create data to save for pheweb
  # clean
  rm(data,data1,data2)
  rm(data_excluded,data_sig)
  rm(pmin,lambda1,lambda,mchisq,mchisq1)

  } else {
    print("Manhattan plot already exists, skip and continue")
  }

  # some cleaning...
  rm(name)
  rm(chr,pheno,pheno_root,subset_name,info_thr,maf_thr,p_col)
  rm(pattern1)
  #

  print(Sys.time())
  print("---------------------------------------------")
  # 
  gc()
}
