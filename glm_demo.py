# this is the demo for how to run the code
import scipy.io as sio
import numpy as np
import os
from scipy.io import savemat
import scipy.stats as stats
from generate_data import generate_SC, generate_FC, generate_coupling
from GLM import generate_covariates, glm

cwd = os.getcwd()
# make your own data dir if not exist
if not os.path.exists(cwd + '/data/'):
    os.makedirs(cwd + '/data/')
data_dir = cwd + '/data/'

if not os.path.exists(cwd + '/glm_res/'):
    	os.makedirs(cwd + '/glm_res/')
glm_dir = cwd + '/glm_res/'

# set the number of the unrelated subjects
# set the number of brain regions in your atlas
unrelated_subj_num = 100
region_num = 392

# generate your own SC, results saved in data dir
generate_SC(unrelated_subj_num, region_num)

# generate your own FC, results saved in data dir
fc_type = 'prec'
generate_FC(unrelated_subj_num, region_num, fc_type)

# generate your SC-FC coupling using the SC and FC, results saved in data dir
generate_coupling(unrelated_subj_num, fc_type)

# GLM analysis, significant value per covariate and confident interval saved in glm dir
generate_covariates(unrelated_subj_num)
glm(unrelated_subj_num, fc_type)



