
import os
import glob
import commands

dat_dir = '/data/clusterfs/lag/projects/lg-ukbiobank/primary_data/imaging_data/'
t1_file = 'T1/T1_unbiased_brain.nii.gz'
subj_list = glob.glob(os.path.join(dat_dir,'*'))
print len(subj_list)
count = 0
for subj in subj_list:
    if(os.path.exists(os.path.join(subj, t1_file))):
        count = count + 1
        sts, out = commands.getstatusoutput('mri_info ' + os.path.join(subj, t1_file))
        out_dat = out.split('\n')
        if out_dat[18].split(': ')[1]!='RAS':
            print subj
        #print out_dat[18].split(': ')[1]


