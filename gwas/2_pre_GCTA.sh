#!/bin/bash

# Compute heritability estimates using GCTA
# #http://cnsgenomics.com/software/gcta/

# to run after GCTA_calUKBv2_imaging.sh (general script to generate GRMs and some test REMLs, also in: GCTA_REML.sh)

# Phenotypes: FS asymmetry phenotypes computed by xiakon
# residualized for: c("assessment_centre","age","sex",paste("PC",1:10,sep=""),"array")

#----------------------------------------------------------------------
# Notes
#----------------------------------------------------------------------
## will need to include further technical and biological confounds, from Elliot et al. bioRxiv: 
## 		technical: scanner x,y,z, position, motion, height, volumetric scaling factor 
## 		biological: height, diastolic and systolic pressure
## others: ICV ???

## but, collider bias? keep in mind
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
resource_dir=/data/workspaces/lag/shared_spaces/Resource_DB/
analysis_dir=/home/xiakon/MPI_workspace/lg-ukbiobank/analysis/xiangzhen/brainTorque/torque_ukb_N40681/GWAS_Dick/
#primary_dir=${analysis_dir}release_v3/cal/QC/
working_dir=${analysis_dir}calc_h2/
UKB_phenos_dir=${analysis_dir}GWAS_out/
working_dir2=/data/clusterfs/lag/users/xiakon/ukb_torque_gwas/
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# run REML - run in multi.q, test running in big.queue. if memory is not so high, just multi.q
#----------------------------------------------------------------------
# copy the data (GRMs) to sge2 directory first, otherwise it won't run in the grid
mkdir -p ${working_dir} ${working_dir}/reml ${working_dir}/logs
mkdir -p ${working_dir2}/grm/ ${working_dir2}/reml ${working_dir2}/pheno ${working_dir2}/logs

# copy genetic data to working_dir2 in sge, if not present yet!
if [ ! -f ${working_dir2}/grm/ukb_cal_snpQC_imagingT1_N33628_clean_rm025_adj.grm.bin ] ; then cp /data/clusterfs/lag/users/xiakon/ukb3/grm/ukb_cal_snpQC_imagingT1_N33628_clean_rm025_adj.grm* ${working_dir2}/grm/ ; fi

#--------------------------------------------
# generate templates, which will be different depending on the value of pheno_root
#--------------------------------------------
echo "
#!/bin/sh
#$ -N REML_multilateral_calUKBv3_Skews
#$ -cwd
#$ -q single.q
#$ -S /bin/bash
#$ -M xiangzhen.kong@mpi.nl
#$ -m as

#----------------------------------------------------------------------
# Define paths and directories
#----------------------------------------------------------------------
working_dir2=/data/clusterfs/lag/users/xiakon/ukb_torque_gwas/
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# get parameters from command line arguments:
root=\$1 # ukb_cal_snpQC_imagingT1_N18057_clean_rm025_adj # Edit
pheno_root=\$2 # 
covs_name=\$3
phenotype=\$4 # PHENOTYPE
#----------------------------------------------------------------------
cd ${working_dir2}
echo 'Phenotype: ' \${phenotype}
#----------------------------------------------------------------------
pheno_file=\${working_dir2}/pheno/\${pheno_root}_residuals4genetics\${covs_name}.table
head -n 1 \$pheno_file | sed 's/ /\n/g' > \${pheno_root}_\${phenotype}.header
# get which column does the phenotype correspond to
pheno_n=\$(grep -nr \$phenotype \${pheno_root}_\${phenotype}.header | awk -F':' '{print \$1}')
# make temporary phenotype file including only : ID1, ID2, phenotype; and no header!
awk -v c1=\$pheno_n '{print \$1,\$2,\$c1}' \$pheno_file | tail -n +2 > \${working_dir2}/pheno/\${pheno_root}_\${phenotype}.pheno
echo 'Check header of phenotype file'
head \${working_dir2}/pheno/\${pheno_root}_\${phenotype}.pheno # just to make sure that it's not the same for all

if [ ! -f \${working_dir2}/reml/\${pheno_root}_reml_\${phenotype}\${covs_name}.hsq ]
then
if [ -f \${working_dir2}/pheno/\${pheno_root}_\${phenotype}.pheno ]
then
echo 'One-grm, quantitative, residuals after adjusting for covariates'
gcta64 --reml --grm \${working_dir2}/grm/\${root} \
 --pheno \${working_dir2}/pheno/\${pheno_root}_\${phenotype}.pheno \
 --thread-num 10 --out \${working_dir2}/reml/\${pheno_root}_reml_\${phenotype}\${covs_name} \
 > \${working_dir2}'logs/'\${pheno_root}'_reml_'\${phenotype}\${covs_name}'.log'
fi
fi

# clean intermediate files
rm \${pheno_root}_\${phenotype}.header 
rm \${working_dir2}/pheno/\${pheno_root}_\${phenotype}.pheno

" > ${analysis_dir}/2_pre_GCTA_Skews_job.sh
#--------------------------------------------

#--------------------------------------------
# Define list of phenotypes: 
#--------------------------------------------

#--------------------------------------------
# Run REML for phenotypes within list 
#--------------------------------------------
# copy phenotype files
cp ${UKB_phenos_dir}/skewsAll_*residuals4genetics*table ${working_dir2}/pheno/
cp ${UKB_phenos_dir}/skewsAll_*residuals4genetics*list ${working_dir2}/pheno/
#cp ${UKB_phenos_dir}/skewscalehand_*residuals4genetics*table ${working_dir2}/pheno/
#cp ${UKB_phenos_dir}/skewscalehand_*residuals4genetics*list ${working_dir2}/pheno/

# run loop, or first just one to test
mkdir -p  ${analysis_dir}/sge_jobs/gcta/
mkdir -p ${working_dir2}/sge_jobs/gcta/
cd ${working_dir2}/sge_jobs/gcta/
dos2unix ${working_dir2}/pheno/skewsAll_*residuals4genetics*table ${working_dir2}/pheno/skewsAll_*residuals4genetics*list
#dos2unix ${working_dir2}/pheno/skewscalehand_*residuals4genetics*table ${working_dir2}/pheno/skewscalehand_*residuals4genetics*list
all_lists=$(ls ${working_dir2}/pheno/skewsAll_*residuals4genetics*list | grep -v NumVert) 

# arguments
#root=ukb_cal_snpQC_imagingT1_N18057_clean_rm025_adj #1
root=ukb_cal_snpQC_imagingT1_N33628_clean_rm025_adj #1
#2: pheno_root, from list file
#3: phenotype, from list file

#test one
#phen=residuals_bankssts
#bash ${analysis_dir}/GCTA_REML_calUKBv2_multilateral.sh $root $pheno_root $phen

#only submit job if output does not exist, neither in working_dir nor working_dir2

#-----------------------------------------------------------------------------
cd ${working_dir2}/sge_jobs/gcta/
# make loop!
for list in ${all_lists} # ${pheno_root}_nonAI_LeftRight_phenotypes.list
 do
 # get covariates name from list name
 pheno_root=$( echo ${list} | awk -F'/pheno/' '{print $2}' | sed 's/.list//g' | awk -F'_residuals4genetics' '{print $1}' ) # CLUSTER_asym_Desikan_MEAN_norm_area_8590_sm00_fsaverage_sym_QC
 covs_name=$(echo ${list} | awk -F'residuals4genetics' '{print $2}' | sed 's/.list//g')

  while read phen
   do
     echo '----------------------------------------'
	 echo ${phen}
     if [ ! -f ${working_dir2}/reml/${pheno_root}_reml_${phen}${covs_name}.hsq ] && [ ! -f ${working_dir}/reml/${pheno_root}_reml_${phen}${covs_name}.hsq ]
     then
      # edit job template for this phenotype
     echo ${root} ${pheno_root} ${covs_name} ${phen}
      ##bash ${analysis_dir}/2_pre_GCTA_Skews_job.sh ${root} ${pheno_root} ${covs_name} ${phen}
      qsub -p -1 ${analysis_dir}/2_pre_GCTA_Skews_job.sh  ${root} ${pheno_root} ${covs_name} ${phen}  # If running with SGE, run first, and run the second time with this line commented.
     else
     echo 'Output for '${phen}' already exists, skip.'
     fi
  done < ${list}
  done


#-----------------------------------------------------------------------------#mkdir -p ${working_dir}
#mkdir -p ${working_dir}
#mkdir -p ${working_dir}reml/ ${working_dir}/logs/
#mv ${working_dir2}reml/skewsAll_* ${working_dir}/reml/
#mv ${working_dir2}logs/skewsAll_* ${working_dir}/logs/
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#cd ${working_dir}/reml/
#hsq_file=$(ls -1 *hsq | grep diff -v)
#paste <(grep "^n" ${hsq_file} ) <(grep "/" ${hsq_file} | awk '{print $2,$3}') <(grep "Pval" ${hsq_file} | awk '{print $2}') > hsq_summary_Skews_v3cal_${covs_name}.table

#-----------------------------------------------------------------------------
## plot estimates: Plot_hsq_results_FSphenos_xiakon.R/3_plot_hsq_results.R

