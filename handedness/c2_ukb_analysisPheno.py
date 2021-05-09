#
import pandas as pd
import numpy as np
import os
import statsmodels.api as sm
from scipy.stats import pearsonr, bartlett, levene

data_file = './ukb_pheno_T1pos.csv'
dat = pd.read_csv(data_file)
print 'Handedness: 1=right, 2=left'
print dat[dat['Handedness']==1]['Handedness'].count()
print dat[dat['Handedness']==2]['Handedness'].count()
print dat[dat['Handedness']==3]['Handedness'].count()

#conf_str = 'plusPCs' # SexAge or plusPCs or plusPos
conf_str = 'plusPos'#'plusPCs' # SexAge or plusPCs or plusPos
#conf_str = 'plusPosScales'#'plusPCs' # SexAge or plusPCs or plusPos
if conf_str == 'SexAge':
    col_covars = ['Sex', 'Age', 'ZAge2'#, 'BrainV'
                 ]
elif conf_str == 'plusPCs':
    col_covars = ['Sex', 'Age','ZAge2', #'BrainV',
                  'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7',
                  'PC8', 'PC9', 'PC10']
elif conf_str == 'plusPos':
    col_covars = ['Sex', 'Age','ZAge2', #'BrainV',
                  'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7',
                  'PC8', 'PC9', 'PC10',
                  'PosX','PosY','PosZ', 'T1SNR', 'T1CNR', 'Center']#+ ['SkewsXY']
else: # plusPosScales
    col_covars = ['Sex', 'Age','ZAge2', 'BrainV',
                  'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7',
                  'PC8', 'PC9', 'PC10',
                  'PosX','PosY','PosZ', 'T1SNR', 'T1CNR', 'Center',
                  'ScalesAvg','ScalesX','ScalesY', 'ScalesZ']

skew_data_file = '../out_2/skewsAll.out'
skew_dat = pd.read_csv(skew_data_file)
#col_skews = ['SkewsXY', 'SkewsXZ'] # 'SkewsYZ'
col_skews = ['ScalesAvg', 'ScalesX', 'ScalesY', 'ScalesZ']   # Do not use plusPosScales for using  this line

# Deal with outliers, 4SD
#TagOutliers = False#True
TagOutliers = True
if TagOutliers:
    thr_xy = skew_dat['SkewsXY'].mean()+4.0*skew_dat['SkewsXY'].std()
    skew_dat.mask(skew_dat['SkewsXY']>thr_xy, inplace=True) # Edit
    skew_dat.mask(skew_dat['SkewsXY']<-thr_xy, inplace=True) # Edit
    thr_xz = skew_dat['SkewsXZ'].mean()+4.0*skew_dat['SkewsXZ'].std()
    skew_dat.mask(skew_dat['SkewsXZ']>thr_xz, inplace=True) # Edit
    skew_dat.mask(skew_dat['SkewsXZ']<-thr_xz, inplace=True) # Edit


## For Handedness ----------------------------------------------------------
hand_dat = pd.merge(left=skew_dat, right=dat, left_on='SID', right_on='SID')
hand_dat.mask(hand_dat['Handedness']==-3, inplace=True)
hand_dat = hand_dat.dropna()
#hand_dat.to_csv('./hand_skew_pheno_T1pos.csv', index=False)

hand_dat['Hand'] = np.nan
hand_dat['Hand'][hand_dat['Handedness']==1] = 1
hand_dat['Hand'][hand_dat['Handedness']==2] = -1
hand_dat['Hand'][hand_dat['Handedness']==3] = 0
#hand_dat.mask(hand_dat['Hand']==0, inplace=True) # Edit
#hand_dat.mask(hand_dat['Hand']==-1, inplace=True) #Edit
#hand_dat.mask(hand_dat['Hand']==1, inplace=True) #Edit
hand_dat = hand_dat.dropna()
hand_dat.to_csv('./hand_skew_pheno_T1pos_clean.csv', index=False)

col_vars = ['Hand']
for col_skew in col_skews:
    print col_skew
    y = hand_dat[col_skew]
    for col_var in col_vars:
        print '-', col_var
        #print pearsonr(hand_dat[col_skew], hand_dat[col_var])
        X = hand_dat[[col_var]+col_covars]
        X = sm.add_constant(X)
        model = sm.OLS(y,X).fit()
        predictions = model.predict(X)
        #print model.summary()
        print 't=', model.tvalues[1], ',p =',model.pvalues[1]

print 'Handedness: 1=right, 2=left, 3=both'
print hand_dat[hand_dat['Handedness']==1]['Handedness'].count()
print hand_dat[hand_dat['Handedness']==2]['Handedness'].count()
print hand_dat[hand_dat['Handedness']==3]['Handedness'].count()

print bartlett(hand_dat[hand_dat['Hand']==-1]['SkewsXY'], hand_dat[hand_dat['Hand']==1]['SkewsXY'])
print bartlett(hand_dat[hand_dat['Hand']==-1]['SkewsXZ'], hand_dat[hand_dat['Hand']==1]['SkewsXZ'])

print levene(hand_dat[hand_dat['Hand']==-1]['SkewsXY'], hand_dat[hand_dat['Hand']==1]['SkewsXY'])
print levene(hand_dat[hand_dat['Hand']==-1]['SkewsXZ'], hand_dat[hand_dat['Hand']==1]['SkewsXZ'])

print hand_dat.groupby('Hand').apply(np.std)

## For Grip Asymmetry --------------------------------------------
dat['Grip'] = dat[['GripL','GripR']].max(axis=1)

z1 = []
z2 = []
z3 = []
z4 = []
for index, row in dat.iterrows():
    if row.Handedness==1:
        z1.append(row.GripR)
        z2.append(row.GripL)
        z3.append(2.0*(row.GripR-row.GripL)/(row.GripR+row.GripL))
        z4.append(2.0*(row.GripR-row.GripL)/(row.GripR+row.GripL))
    elif row.Handedness==2:
        z1.append(row.GripL)
        z2.append(row.GripR)
        z3.append(2.0*(row.GripL-row.GripR)/(row.GripR+row.GripL))
        z4.append(2.0*(row.GripR-row.GripL)/(row.GripR+row.GripL))
    else:
        z1.append(np.nan)
        z2.append(np.nan)
        z3.append(np.nan)
        z4.append(np.nan)
dat['yGrip'] = z1 # dominance hand
dat['nGrip'] = z2
dat['GripAsyAbs'] = z3
dat['GripAsy'] = z4 # for the main analysis, grip strength asymmetry !!!

print dat.shape

col_vars = ['GripAsy']
print 'GripAsy #: ', dat['GripAsy'].count()

grip_dat = pd.merge(left=skew_dat, right=dat, left_on='SID', right_on='SID')
grip_dat = grip_dat.dropna()
grip_dat.to_csv('./grip_skew_pheno_T1pos.csv', index=False)

col_vars = ['GripAsy']
for col_skew in col_skews:
    print col_skew
    y = grip_dat[col_skew]
    for col_var in col_vars:
        print '-', col_var
        #print pearsonr(grip_dat[col_skew], grip_dat[col_var])
        X = grip_dat[[col_var]+col_covars]
        X = sm.add_constant(X)
        model = sm.OLS(y,X).fit()
        predictions = model.predict(X)
        #print model.summary()
        print 't=', model.tvalues[1], ',p =',model.pvalues[1]



"""
# for subject info of interest
info_fields = {'SID': 'eid', 
		'Sex': '31-0.0', 
		'Age': '21003-2.0', 
		'Handedness': '1707-2.0', 
		'GripL':'46-2.0',
		'GripR':'47-2.0',
		'Height':'12144-2.0',
		'Weight':'12143-2.0'}
# N = 9932 for subcortical volume, overall white/gray matter
# N = 9931 for cortical volume
# N = 8839 for DTI measures  !!!!!!!!!!!!
# N = 9306 for rfMRI matrix
#
# N = 11671 for handedness
# N = 12551 for age for imaging
# N = 502629 for sex for all
# here we focus on subject info with DTI meaures, take collumn 25062-2.0, ...
# ... Mean FA in corticospinal tract on FA skeleton (right), as the selection rule
# here we focus on subject info with T1 data, take collumn 25006-2.0, ...
my_index = dat['25006-2.0'].notnull()
col_names = ['SID', 'Sex', 'Age', 'Handedness', 'GripL', 'GripR', 'Height', 'Weight']
col_list = [info_fields['SID'], info_fields['Sex'], info_fields['Age'], 
		info_fields['Handedness'],
		info_fields['GripL'], info_fields['GripR'],
		info_fields['Height'], info_fields['Weight']]
info_dat = dat[col_list][my_index]
info_dat.columns = col_names
info_dat.to_csv('../UKB_Dat/ukb_pheno.csv', index=False)
"""
