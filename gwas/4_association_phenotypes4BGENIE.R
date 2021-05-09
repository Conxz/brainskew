# association_phenotypes4BGENIE.R
# date: 13.10.2017; edited: 09.05.2018; edited 23.05.2018; edited 26.09.2018 to adapt for FSphenos_xiakon

# create phenotype input files for BGENIE (which do not have sample IDs, and need to match the order within the BGEN file)
# genetic data: v3; imputed data, 
## includes chromosome X, and XY, but sample files are different for these chromosomes, so new input files will be required


#------------------------------------
# Get command line arguments
args = commandArgs(trailingOnly=TRUE)

# args=c("/data/clusterfs/lag/users/amacar/ukb/input/imp/ukb_imp_chr19_v3_imagingT1_N12245.sample",
#        "A",
#        "GWAS",
#        "imagingT1_N12245")

# get sample file as argument
sample_file=args[1]
chr_type=args[2] # A for autosomes, X or XY
pheno_root=args[3] # ukb phenotype batch names, or type of analysis
subset_name=args[4] # subset of data that will be included,, in case not samples are in; should reflect sample_file info

# FS_type=args[5]
#------------------------------------


options(stringsAsFactors = FALSE)
#------------------------------------
# Set directories, dependent on  system
if (Sys.info()['sysname']=='Windows') {dir="P://workspaces/"} else {dir="/data/workspaces/lag/workspaces/"}
imaging_dir=paste(dir,"lg-ukbiobank/analysis/xiangzhen/brainTorque/torque_ukb_N40681/",pheno_root,"/",sep="")
out_dir_base=paste(dir,"lg-ukbiobank/analysis/xiangzhen/brainTorque/torque_ukb_N40681/GWAS_Dick/GWAS_out2/",sep="")
out_dir=paste(out_dir_base,"imp_v3_bgenie/",sep="")
dir.create(file.path(out_dir_base,"imp_v3_bgenie"),showWarnings = FALSE)

#------------------------------------
# Sample files,  read
#------------------------------------

## read
sample<-read.table(sample_file,header=TRUE) # actually, use v3: they have the same number of individuals for the autosomes
# some checks
table(sample$missing)
#table(sample$sex)


# remove first line which is only 0
sample[1,]
sample<-sample[-1,1:2]

#------------------------------------
# Phenotype file for imaging phenotypes, output from: 
## ukb9246_phenos_multilateral.R
#------------------------------------
out_file=paste(out_dir,pheno_root,"_sample_",subset_name,"_phenos4BGENIE.table",sep="")
# Run all phenos at once
files<-list.files(out_dir_base,pattern="*_noBioCovs_ZAge.table")

for (f in files){
  measure<-sapply(strsplit(f,"_"),"[[",1)
  asy<-sapply(strsplit(f,"_"),"[[",3)
  #pheno_imaging_file<-paste(imaging_dir,f,sep="")
  pheno_imaging_file<-paste(out_dir_base,f,sep="")
  #print(pheno_imaging_file)
  #  read data
  pheno_imaging<-read.table(pheno_imaging_file,header=TRUE)
  # get phenos:
  phenos<-colnames(pheno_imaging)[grep("residuals",colnames(pheno_imaging))]
  # subset and rename
  colnames(pheno_imaging)<-gsub("Asy",paste(measure,"_Asy",sep=""),colnames(pheno_imaging))
  colnames(pheno_imaging)<-gsub("Asy",asy,colnames(pheno_imaging)) # add info about Asy or HealthyAsy
  # combine all
  if (f==files[1]){
    all_imaging<-pheno_imaging
  } else {
      all_imaging<-merge(all_imaging,pheno_imaging,stringAsFactors=FALSE,all=TRUE)
      }
  # clean
  rm(measure,asy,pheno_imaging_file,pheno_imaging,phenos)
}
rm(files,f)

# select phenos
phenos<-colnames(all_imaging)[grep("residuals",colnames(all_imaging))]
#phenos_select<-phenos[grep("Skews",phenos)]
phenos_select<-phenos[grep("residuals",phenos)]

# combine pheno files
if (length(grep("X",chr_type))>0) {
  out_file2<-gsub("sample",paste("sample",chr_type,sep=""),out_file)
} else {
  out_file2<-out_file 
  }

if (file.exists(out_file2)==FALSE){
  # make sure that the pheno file matches the order from the sample
  pheno_samples<-merge(sample,all_imaging,by.x=c("ID_1"),by.y=c("ID1"),all.x=TRUE,stringsAsFactors=FALSE) #,
  pheno_samples$array<-as.character(pheno_samples$array)
  pheno_samples$batch<-as.character(pheno_samples$batch)
  # check if order is the same
  table(sample$ID_1==pheno_samples$ID_1)
  pheno_samples<-pheno_samples[match(sample$ID_1,pheno_samples$ID_1),]
  table(sample$ID_1==pheno_samples$ID_1)
  
  # missing data is NA, need to specify this when running BGENIE, replace to -999
  pheno_samples[is.na(pheno_samples)]<-(-999)
  # or alternatively use the --miss parameter to code for missingness (NA)
  
  # save file with phenotypes and covariates
  # include just one ID, to be able to double check
  write.table(pheno_samples[,c("ID_1",phenos)],file=out_file2,row.names=FALSE,quote=FALSE)
  
  rm(pheno_samples)
}

