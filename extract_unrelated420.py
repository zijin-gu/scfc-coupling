import scipy.io as sio
import numpy as np
import os

cwd = os.getcwd()
data_dir = cwd + '/data/'

all_subj  = np.genfromtxt(data_dir + 'subjects_all997.txt', dtype = 'int')
unrelated_subj = np.genfromtxt(data_dir + 'subjects_unrelated420.txt', dtype = 'int')

roi_number = 392

def extract_sc(sc_type):
	
	sc = np.zeros([unrelated_subj.shape[0],roi_number,roi_number])

	sc_all = sio.loadmat(data_dir + 'sc_ifod2act_cc400_997subj.mat')

	index = all_subj.searchsorted(unrelated_subj)
	i = 0
	for idx in index:
		sc[i] = sc_all[sc_type][0][idx]
		i += 1
		
	mdic = {"C": sc}
	sio.savemat(data_dir + 'cc400_SC_' + sc_type + '_unrelated420.mat', mdic)

def extract_fc_corr(filter_type):

	fc_dir = data_dir + 'fmriclean_concat_' + filter_type + '_cc400_gsr_FCcov_997subj/'

	fc_corr = np.zeros([unrelated_subj.shape[0],roi_number,roi_number])
	
	n = 0
	for subj in unrelated_subj:
	    fccov = sio.loadmat(fc_dir + str(subj) + '_concat_fmriclean_' + filter_type + '_cc400_gsr_FCcov.mat')
	    fc_corr[n] = fccov['C']
	    n += 1

	mdic = {"C": fc_corr}
	sio.savemat(data_dir + 'cc400_FCcorr_concat_' + filter_type + '_unrelated420.mat', mdic)

def extract_fc_prec(filter_type):

	fc_dir = data_dir + 'fmriclean_concat_' + filter_type + '_cc400_gsr_FCcov_997subj/'
	fc_corr = np.zeros([unrelated_subj.shape[0],roi_number,roi_number])

	n = 0
	for subj in unrelated_subj:
	    fccov = sio.loadmat(fc_dir + str(subj) + '_concat_fmriclean_' + filter_type + '_cc400_gsr_FCcov.mat')
	    fc_corr[n] = fccov['C']
	    n += 1

	# regularized prec 
	min_rmse = 1e5
	opt_gamma = 0
	for gamma in np.linspace(0,1,100):
	    inverse = []
	    reg_inv = []
	    for i in range(unrelated_subj.shape[0]):
	        np.fill_diagonal(fc_corr[i],1)
	        inverse.append(np.linalg.inv(fc_corr[i]))
	        reg_inv.append(np.linalg.inv(fc_corr[i] + gamma*np.eye(roi_number)))
	    group_prec = np.mean(inverse, axis=0)
	    diff = []
	    for i in range(unrelated_subj.shape[0]):
	        diff.append(np.linalg.norm(reg_inv[i][np.triu_indices(roi_number,1)] - group_prec[np.triu_indices(roi_number,1)]))
	    rmse = np.mean(diff)
	    if rmse < min_rmse:
	        min_rmse = rmse
	        opt_gamma = gamma

	fc_prec = np.zeros([unrelated_subj.shape[0],roi_number,roi_number])
	for i in range(unrelated_subj.shape[0]):
	    fc_prec[i] = np.linalg.inv(fc_corr[i] + opt_gamma*np.eye(roi_number))

	mdic = {"C": fc_prec, "gamma": opt_gamma}
	sio.savemat(data_dir + 'cc400_FCprec_concat_' + filter_type + '_unrelated420.mat', mdic)
