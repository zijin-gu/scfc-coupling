# load data

import scipy.io as sio
import numpy as np
from scipy.io import savemat

scfc_dir = './data/scfc_data/'
alldata_dir = './data/allsubjects/'
unrelated420_dir = './data/unrelated420/'

with open(alldata_dir + 'subjects_rfMRI_dMRI_complete_997.txt') as all_subj:
	    all_subjects = all_subj.readlines()

with open(unrelated420_dir + 'subjects_unrelated420_scfc.txt') as unrelated_subj:
    unrelated_subjects = unrelated_subj.readlines()

all_subjects = [ int(x) for x in all_subjects ]
unrelated_subjects = [ int(x) for x in unrelated_subjects ]

def extract_sc():
	
	sc_all = sio.loadmat(sc_dir + 'sc_ifod2act_cc400_997subj.mat')

	sc_orig = np.zeros([len(unrelated_subjects),392,392])
	sc_sift2 = np.zeros([len(unrelated_subjects),392,392])
	sc_volnorm = np.zeros([len(unrelated_subjects),392,392])
	sc_sift2volnorm = np.zeros([len(unrelated_subjects),392,392])

	# remove the related subjects
	for n in range(len(unrelated_subjects)):
	    for m in range(len(all_subjects)):
	        if all_subjects[m] == unrelated_subjects[n]:
	            sc_orig[n] = sc_all['orig'][0][m]
	            sc_sift2[n] = sc_all['sift2'][0][m]
	            sc_volnorm[n] = sc_all['volnorm'][0][m]
	            sc_sift2volnorm[n] = sc_all['sift2volnorm'][0][m]
	            break    
	
	mdic = {"C": sc_sift2volnorm}
	savemat(unrelated420_dir + "cc400_SC_sift2volnorm_unrelated420.mat", mdic)


def extract_FCcorr_nofit():
	fcgsr_all = sio.loadmat(scfc_dir + 'fc_gsr_cc400_997subj.mat')
	fcgsr = np.zeros([len(unrelated_subjects),392,392])

	# remove the related subjects
	for n in range(len(unrelated_subjects)):
	    for m in range(len(all_subjects)):
	        if all_subjects[m] == unrelated_subjects[n]:
	            fcgsr[n] = fcgsr_all['fc'][0][m]
	            break
	
	mdic = {"C": fcgsr}
	savemat(unrelated420_dir + "cc400_FCcorr_nofit_unrelated420.mat", mdic)


def extract_FCprec_nofit():
	FCgsr_prec_unrelated420 = np.zeros([len(unrelated_subjects),392,392])
	fcdir = scfc_dir + 'fmriclean_nofilt_cc400_gsr_FCprec_997subj/'
	
	n = 0
	for subj in unrelated_subjects:
	    fcprec_gsr1LR = sio.loadmat(fcdir + str(subj) + '_rfMRI_REST1_LR_fmriclean_nofilt_cc400_gsr_FCprec.mat')
	    fcprec_gsr1RL = sio.loadmat(fcdir + str(subj) + '_rfMRI_REST1_RL_fmriclean_nofilt_cc400_gsr_FCprec.mat')
	    fcprec_gsr2LR = sio.loadmat(fcdir + str(subj) + '_rfMRI_REST2_LR_fmriclean_nofilt_cc400_gsr_FCprec.mat')
	    fcprec_gsr2RL = sio.loadmat(fcdir + str(subj) + '_rfMRI_REST2_RL_fmriclean_nofilt_cc400_gsr_FCprec.mat')
	    fcprec = (fcprec_gsr1LR['C'] + fcprec_gsr1RL['C'] + fcprec_gsr2LR['C'] + fcprec_gsr2RL['C'])/4
	    FCgsr_prec_unrelated420[n] = fcprec
	    n += 1
	
	mdic = {"C": FCgsr_prec_unrelated420}
	savemat(unrelated420_dir + "cc400_FCprec_nofit_unrelated420.mat", mdic)

# filter_type = nofilt, bpf, hpf
# extract the corr and prec version for 420 unrelated subjects
def extract_concat(filter_type):
	FCgsr_cov_unrelated420 = np.zeros([len(unrelated_subjects),392,392])
	fcdir = scfc_dir + 'fmriclean_concat_hpf_cc400_gsr_FCcov_997subj/'
	n = 0
	for subj in unrelated_subjects:
	    fccov = sio.loadmat(fcdir + str(subj) + '_concat_fmriclean_' + filter_type + '_cc400_gsr_FCcov.mat')
	    FCgsr_cov_unrelated420[n] = fccov['C']
	    n += 1

	mdic = {"C": FCgsr_cov_unrelated420}
	savemat(unrelated420_dir + 'cc400_FCcorr_concat_' + filter_type + '_unrelated420.mat', mdic)

	# regularized prec 
	min_rmse = 1e5
	opt_gamma = 0
	for gamma in np.linspace(0,1,100):
	    inverse = []
	    reg_inv = []
	    for i in range(len(unrelated_subjects)):
	        np.fill_diagonal(FCgsr_cov_unrelated420[i],1)
	        inverse.append(np.linalg.inv(FCgsr_cov_unrelated420[i]))
	        reg_inv.append(np.linalg.inv(FCgsr_cov_unrelated420[i] + gamma*np.eye(392)))
	    group_prec = np.mean(inverse,axis=0)
	    diff = []
	    for i in range(len(unrelated_subjects)):
	        diff.append(np.linalg.norm(reg_inv[i][np.triu_indices(392,1)] - group_prec[np.triu_indices(392,1)]))
	    rmse = np.mean(diff)
	    if rmse < min_rmse:
	        min_rmse = rmse
	        opt_gamma = gamma

	FCgsr_prec = np.zeros([len(unrelated_subjects),392,392])
	for i in range(len(unrelated_subjects)):
	    FCgsr_prec[i] = np.linalg.inv(FCgsr_cov_unrelated420[i] + opt_gamma*np.eye(392))

	mdic = {"C": FCgsr_prec}
	savemat(unrelated420_dir + 'cc400_FCprec_concat_' + filter_type + '_unrelated420.mat', mdic)
