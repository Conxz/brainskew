
import os
from glob import glob
import numpy as np

#datDir = '/data/clusterfs/lag/users/xiakon/vol_asym/dat/'
datDir = '/data/clusterfs/lag/users/xiakon/torque/dat'

sidListFile = '../doc/sidList.txt'
sidList = [sid.strip() for sid in open(sidListFile)]
print len(sidList), 'subjects'

fileName = 'T1/transforms/T1_to_MNIsymm_linear_dof12.mat'

outDir = '../out_1'
if not os.path.exists(outDir):
    os.mkdir(outDir)

#docDir = 'doc'
#if not os.path.exists(docDir):
#    os.mkdir(docDir)
#sidAll = []
#sidNull = []

#run_sge_str = 'fsl_sub -q single.q -l ./sgelog ' # seems not working well!!!

for sid in sidList:
    sidPath = os.path.join(datDir, sid)
    matFile = os.path.join(sidPath, fileName)
    #sid = os.path.basename(sidPath)
    #print sid
    if os.path.exists(matFile):
        outFile = os.path.join(outDir, sid+'_avscale.out')
        if not os.path.exists(outFile):
            #print matFile
            print 'avscale ' + matFile + ' > ' + outFile
            #os.system(run_sge_str + 'avscale ' + matFile + ' > ' + outFile) # seems not working well!!!
            os.system('avscale ' + matFile + ' > ' + outFile) # quick enough.
            #sidAll.append(sid)
    else:
        #sidNull.append(sid)
        print sidPath

#np.savetxt(os.path.join(docDir, 'sidList.txt'), sidAll, fmt='%s')
#np.savetxt(os.path.join(docDir, 'sidNull.txt'), sidNull, fmt='%s')

