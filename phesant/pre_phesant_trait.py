
import pandas as pd
import os

csv_file = '../out_2/skewsAll.out'
csv_dat = pd.read_csv(csv_file)

trait_list = ['SkewsXY','SkewsXZ','SkewsYZ']

out_dir = './'

for trait in trait_list:
    print trait
    out_file = os.path.join(out_dir, trait+'_PHESANT.csv')
    csv_dat[['SID']+[trait]].to_csv(out_file, index=False)

