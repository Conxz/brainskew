
import os
from glob import glob
import numpy as np

sidListFile = '../doc/sidList.txt'
sidList = [sid.strip() for sid in open(sidListFile)]

niiDir = '/data/clusterfs/lag/users/xiakon/torque/dat'
niiFileName = 'T1/T1_brain_to_MNIsymm_dof9.nii.gz'

outDir = '/data/clusterfs/lag/users/xiakon/torque/dat'
z_slice = 73
z_value = str(round(z_slice/192.0, 2))
figFileName_z = 'dof9_z'+str(z_slice)+'.png'

y_slice = 105
y_value = str(round(y_slice/228.0, 2))
figFileName_y = 'dof9_y'+str(y_slice)+'.png'

x_slice = 90
x_value = str(round(x_slice/192.0, 2))
figFileName_x = 'dof9_x'+str(x_slice)+'.png'

docDir = '../doc'

run_sge_str = 'fsl_sub -q single.q -l ./sgelog '
#run_sge_str = ''
sidError = []
for sid in sidList:
    print sid
    niiFile = os.path.join(niiDir, sid, niiFileName)
    
    if os.path.exists(niiFile):
        sidDir = os.path.join(outDir, sid)
        if not os.path.exists(sidDir):
            os.mkdir(sidDir)
        outFileDir = os.path.join(sidDir, 'T1')
        if not os.path.exists(outFileDir):
            os.mkdir(outFileDir)
        
        outFile_z = os.path.join(outFileDir, figFileName_z)
        outFile_y = os.path.join(outFileDir, figFileName_y)
        outFile_x = os.path.join(outFileDir, figFileName_x)
        
        if not os.path.exists(outFile_z):
            os.system(run_sge_str + 'slicer ' + niiFile + ' -u -z ' + str(z_value) + ' ' + outFile_z)
        if not os.path.exists(outFile_y):
            os.system(run_sge_str + 'slicer ' + niiFile + ' -u -y ' + str(y_value) + ' ' + outFile_y)
        if not os.path.exists(outFile_x):
            os.system(run_sge_str + 'slicer ' + niiFile + ' -u -x ' + str(x_value) + ' ' + outFile_x)
        #print run_sge_str + 'slicer ' + niiFile + ' -u -z ' + str(z_value) + ' ' + outFile
    else:
        sidError.append(sid)

np.savetxt(os.path.join(docDir, 'sidError_dof9_slicer.txt'), sidError, fmt='%s')

