
# for running PHESANT 
import os

phesant_dir = '/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/xiangzhen/ukb/PHESANT'
phesant_r = 'WAS/phenomeScan.r'
phesant_script = os.path.join(phesant_dir, phesant_r)

variablelist_file = os.path.join(phesant_dir, 'variable-info/outcome-info.tsv')
datacoding_file = os.path.join(phesant_dir, 'variable-info/data-coding-ordinal-info.txt')

project_dir = '/data/workspaces/lag/workspaces/lg-ukbiobank/analysis/xiangzhen/brainTorque/torque_ukb_N40681/PHESANT'
pheno_file = os.path.join(project_dir, 'ukb_pheno_dat4skew.csv') # edit
user_id = 'SID'

conf_file = os.path.join(project_dir, 'confounders_scannerPosition_t1.csv') # edit

trait_dir = project_dir
#trait_list = ['SkewsXY', 'SkewsXZ', 'SkewsYZ']
trait_list = ['SkewsXY', 'SkewsXZ']
#trait_list = ['SkewsYZ']

tmp_dir = os.getcwd()
os.chdir(os.path.join(phesant_dir, 'WAS'))
for trait in trait_list:
    print '----------------------------------'
    print '--------------'+trait+'----------------'
    print '---->Running PHESANT ....'
    os.chdir(os.path.join(phesant_dir, 'WAS'))
    out_dir = os.path.join(trait_dir, 'res_'+trait+'_conf_scanT1/')#'_conf_scanT1_gs/') # gs for genetic and sensitivity # Note that the last '/' is a must
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    traitofinterest = trait
    trait_file = os.path.join(trait_dir, trait+'_PHESANT.csv')
    print pheno_file, trait_file, variablelist_file, datacoding_file, traitofinterest, conf_file
    cmd_str = phesant_script+' --phenofile='+pheno_file+' --traitofinterestfile='+trait_file+' --variablelistfile='+\
              variablelist_file+' --datacodingfile='+datacoding_file+' --traitofinterest='+traitofinterest+\
              ' --confounderfile='+conf_file +\
              ' --genetic=TRUE --sensitivity' +\
              ' --resDir='+out_dir+' --userId='+user_id
    cmd_str = 'Rscript ' + cmd_str
    print cmd_str
    os.system(cmd_str)

    print '---->Combining results ...'
    os.chdir(os.path.join(phesant_dir, 'resultsProcessing'))
    cmd_str2 = 'mainCombineResults.r --resDir='+out_dir + ' --variablelistfile='+variablelist_file
    cmd_str2 = 'Rscript ' + cmd_str2
    print cmd_str2
    os.system(cmd_str2)

    print '---->Visualizing results ....'
    os.chdir(os.path.join(phesant_dir, 'PHESANT-viz/bin'))
    web_dir = os.path.join(out_dir, 'web')
    if not os.path.exists(web_dir):
        os.mkdir(web_dir)
    if not os.path.exists(os.path.join(web_dir, trait+'.html')):
        web_file_src = os.path.join(phesant_dir, 'PHESANT-viz/web/force-collapsible-real.html')
        cp_cmd = 'cp ' + web_file_src + ' ' + os.path.join(web_dir, trait+'.html')
        os.system(cp_cmd)
    pheno_res_file = os.path.join(out_dir, 'results-combined.txt')
    json_file = os.path.join(web_dir, 'java-json.json')
    cmd_str3 = 'java -cp .:../jar/json-simple-1.1\ 2.jar ResultsToJSON \"'+ pheno_res_file + '\" \"../node-positions.csv\" \"'+json_file+'\"'
    print cmd_str3
    os.system(cmd_str3)
    
os.chdir(tmp_dir)

print '--------------Finished----------------'
print '----------------------------------'
