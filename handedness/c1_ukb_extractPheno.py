#
import pandas as pd
import numpy as np
import os

dat_dir = '../PHESANT' # data file path
dat_file = os.path.join(dat_dir, 'ukb_pheno_dat4skew.csv')
dat = pd.read_csv(dat_file)

# for subject info of interest
info_fields = {'SID': 'SID', 
		'Sex': 'x31_0_0', 
		'Center': 'x54_2_0', 
		'Age': 'x21003_2_0', 
		#'Handedness': 'x1707_0_0',#'x1707_2_0', 
		'Handedness0': 'x1707_0_0',#'x1707_2_0', 
		'Handedness1': 'x1707_1_0',#'x1707_2_0', 
		'Handedness2': 'x1707_2_0',#'x1707_2_0', 
		'GripL':'x46_2_0',
		'GripR':'x47_2_0',
		'Height':'x12144_2_0',
		'Weight':'x12143_2_0',
		'PC1':'x22009_0_1',
		'PC2':'x22009_0_2',
		'PC3':'x22009_0_3',
		'PC4':'x22009_0_4',
		'PC5':'x22009_0_5',
		'PC6':'x22009_0_6',
		'PC7':'x22009_0_7',
		'PC8':'x22009_0_8',
		'PC9':'x22009_0_9',
		'PC10':'x22009_0_10', 
                'PosX':'x25756_2_0',
                'PosY':'x25757_2_0',
                'PosZ':'x25758_2_0',
                'T1CNR':'x25734_2_0',
                'T1SNR':'x25735_2_0',
                'BrainV':'x25010_2_0'}

# here we focus on subject info with T1 data, take collumn 25006-2.0, ...
my_index = dat['x25006_2_0'].notnull()
#col_names = ['SID', 'Center', 'Sex', 'Age', 'Handedness', 'GripL', 'GripR', 'Height', 'Weight',
col_names = ['SID', 'Center', 'Sex', 'Age', 'Handedness0', 'Handedness1','Handedness2','GripL', 'GripR', 'Height', 'Weight',
             'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9','PC10', 
             'PosX', 'PosY', 'PosZ', 'T1CNR','T1SNR', 'BrainV']
col_list = [info_fields['SID'], info_fields['Center'],
                info_fields['Sex'], info_fields['Age'], 
		#info_fields['Handedness'],
		info_fields['Handedness0'],info_fields['Handedness1'],info_fields['Handedness2'],
		info_fields['GripL'], info_fields['GripR'],
		info_fields['Height'], info_fields['Weight'],
		info_fields['PC1'], info_fields['PC2'],
		info_fields['PC3'], info_fields['PC4'],
		info_fields['PC5'], info_fields['PC6'],
		info_fields['PC7'], info_fields['PC8'],
		info_fields['PC9'], info_fields['PC10'],
		info_fields['PosX'], info_fields['PosY'],
		info_fields['PosZ'], 
                info_fields['T1CNR'], info_fields['T1SNR'],
                info_fields['BrainV']
            ]
info_dat = dat[col_list][my_index]
info_dat.columns = col_names

info_dat['ZAge2'] = ((info_dat['Age'] - info_dat['Age'].mean())/info_dat['Age'].std())**2
info_dat['Handedness'] = info_dat['Handedness2'].combine_first(info_dat['Handedness1']).combine_first(info_dat['Handedness0'])
info_dat = info_dat.drop(['Handedness0', 'Handedness1', 'Handedness2'], axis=1)

info_dat.to_csv('./ukb_pheno_T1pos.csv', index=False)

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
