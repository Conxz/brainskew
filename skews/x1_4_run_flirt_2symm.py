
import os
from glob import glob
import numpy as np

datDir = '/data/clusterfs/lag/projects/lg-ukbiobank/primary_data/imaging_data/'

#sidListFile = './doc/sidList.txt'
#sidList = [sid.strip() for sid in open(sidListFile)]
sidList = glob(os.path.join(datDir, '*'))
sidList = [os.path.basename(sidPath) for sidPath in sidList]

#srcFileName = 'T1/T1_brain.nii.gz'
srcFileName = 'T1/T1_unbiased_brain.nii.gz' # Corrected. use data after bias filed correction.
refFile = '../doc/mni_icbm152_t1_tal_nlin_sym_09c_brain.nii.gz'

outDir = '/data/clusterfs/lag/users/xiakon/torque/dat'
if not os.path.exists(outDir):
    os.mkdir(outDir)

outFileName0 = 'T1_brain_to_MNIsymm_dof6.nii.gz'
outFileName = 'T1_brain_to_MNIsymm_dof9.nii.gz'
outFileName2 = 'T1_brain_to_MNIsymm_dof12.nii.gz'
outMatName0 = 'T1_to_MNIsymm_linear_dof6.mat'
outMatName = 'T1_to_MNIsymm_linear_dof9.mat'
outMatName2 = 'T1_to_MNIsymm_linear_dof12.mat'

docDir = '../doc'

run_sge_str = 'fsl_sub -q single.q -l ./sgelog '
sidError = []
sidAll = []

for sid in sidList:
    print sid
    srcFile = os.path.join(datDir, sid, srcFileName)
    
    if os.path.exists(srcFile):
        sidAll.append(sid)
        sidDir = os.path.join(outDir, sid)
        if not os.path.exists(sidDir):
            os.mkdir(sidDir)
        outFileDir = os.path.join(sidDir, 'T1')
        if not os.path.exists(outFileDir):
            os.mkdir(outFileDir)
        outMatDir = os.path.join(outFileDir, 'transforms')
        if not os.path.exists(outMatDir):
            os.mkdir(outMatDir)
        outFile0 = os.path.join(outFileDir, outFileName0)
        outFile = os.path.join(outFileDir, outFileName)
        outFile2 = os.path.join(outFileDir, outFileName2)
        
        outMat0 = os.path.join(outMatDir, outMatName0)
        outMat = os.path.join(outMatDir, outMatName)
        outMat2 = os.path.join(outMatDir, outMatName2)
        
        if not os.path.exists(outMat0): # note to edit the outMat file name
            os.system(run_sge_str + 'flirt -dof 6 -in ' + srcFile + ' -ref ' + refFile + ' -out ' + outFile0 + ' -omat ' + outMat0) 
            #os.system(run_sge_str + 'flirt -dof 9 -in ' + srcFile + ' -ref ' + refFile + ' -out ' + outFile + ' -omat ' + outMat) 
            #os.system(run_sge_str + 'flirt -dof 12 -in ' + srcFile + ' -ref ' + refFile + ' -out ' + outFile2 + ' -omat ' + outMat2) # Main
    else:
        sidError.append(sid)

#np.savetxt(os.path.join(docDir, 'sidError_flirtdof9_12.txt'), sidError, fmt='%s')
#np.savetxt(os.path.join(docDir, 'sidList.txt'), sidAll, fmt='%s')

