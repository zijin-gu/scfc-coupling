from extract_unrelated420 import extract_sc, extract_fc_corr, extract_fc_prec
from gen_coupling import gen_scfc_coupling
from GLM import glm
from extract_herit941 import extract_herit_sc, extract_herit_fc, gen_herit_scfc_couping, standardize
import scipy.io as sio
import numpy as np
import os
import scipy.stats as stats

cwd = os.getcwd()
data_dir = cwd + '/data/'

all_subj = np.genfromtxt(data_dir + 'subjects_all997.txt', dtype='int')
unrelated_subj = np.genfromtxt(
    data_dir + 'subjects_unrelated420.txt', dtype='int')

roi_number = 392
opt_gamma = 0.3

def extract_data():
    extract_sc('sift2volnorm')
    extract_fc_corr('hpf')
    extract_fc_prec('hpf')


extract_data()
gen_scfc_coupling('sift2volnorm', 'prec', 'hpf')
gen_scfc_coupling('sift2volnorm', 'corr', 'hpf')
glm('prec', 'hpf')
glm('corr', 'hpf')
extract_herit_sc('sift2volnorm')
extract_herit_fc(opt_gamma)
standardize('sc')
for fc_type in ['prec', 'corr']:
	for sess in range(1,5):
		gen_herit_scfc_couping('sift2volnorm', fc_type, sess)
		standardize('fc', fc_type, sess)
		standardize('cp', fc_type, sess)

