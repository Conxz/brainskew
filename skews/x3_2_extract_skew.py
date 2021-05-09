
from glob import glob
import os
import numpy as np
import pandas as pd

outSrc = '../out_1'
outList = glob(os.path.join(outSrc, '*.out'))

outDir = '../out_2'
if not os.path.exists(outDir):
    os.mkdir(outDir)
outFile = os.path.join(outDir, 'skewsAll.out')

sidList = []
skewsXY = []
skewsXZ = []
skewsYZ = []
scalesX = []
scalesY = []
scalesZ = []
avgScale = []

for srcFile in outList:
    sid = os.path.basename(srcFile).split('_')[0]
    sidList.append(sid)

    f = open(srcFile)
    outDat = f.read()
    f.close()
    lineSkews = outDat.split('\n')[8]
    lineScales = outDat.split('\n')[6]
    lineAvgScales = outDat.split('\n')[10]
    
    # for skews
    strSkews = lineSkews.split('=')[0].strip()
    datSkews = lineSkews.split('=')[1].strip().split(' ')
    #print strSkews, datSkews
    skewsXY.append(datSkews[0])
    skewsXZ.append(datSkews[1])
    skewsYZ.append(datSkews[2])
    
    # for scales
    strScales = lineScales.split('=')[0].strip()
    datScales = lineScales.split('=')[1].strip().split(' ')
    #print strScales, datScales
    scalesX.append(datScales[0])
    scalesY.append(datScales[1])
    scalesZ.append(datScales[2])
    
    # for average scaling
    strAvgScales = lineAvgScales.split('=')[0].strip()
    datAvgScales = lineAvgScales.split('=')[1].strip()
    #print strAvgScales, datAvgScales
    avgScale.append(datAvgScales)

outDF = pd.DataFrame({'SID': sidList,
                      'SkewsXY': skewsXY, 
                      'SkewsXZ': skewsXZ,
                      'SkewsYZ': skewsYZ, 
                      'ScalesX': scalesX,
                      'ScalesY': scalesY, 
                      'ScalesZ': scalesZ,
                      'ScalesAvg': avgScale})
outDF.to_csv(outFile, index=False, sep=',')

