import scipy.io as sio
import numpy as np
import os
from scipy.io import savemat, loadmat
import scipy.stats as stats

cwd = os.getcwd()
# make your own data dir if not exist
if not os.path.exists(cwd + '/data/'):
	os.makedirs(cwd + '/data/')
data_dir = cwd + '/data/'

def generate_SC(unrelated_subj_num, region_num):

	sc = np.random.randn(unrelated_subj_num, region_num, region_num)

	mdic = {'C': sc}
	savemat(data_dir + 'SC_unrelated' + str(unrelated_subj_num) + '.mat', mdic)

def generate_FC(unrelated_subj_num, region_num, fc_type):
	# fc_type has to be either "corr" or "prec"

	if fc_type == 'corr':
		fc_corr = np.random.randn(unrelated_subj_num, region_num, region_num)

		mdic = {'C': fc_corr}
		savemat(data_dir + 'FC_' + fc_type + '_unrelated' + str(unrelated_subj_num) + '.mat', mdic)

	elif fc_type == 'prec':

		fc_corr = np.random.randn(unrelated_subj_num, region_num, region_num)

		min_rmse = 1e5
		opt_gamma = 0
		for gamma in np.linspace(0,1,100):
			inverse = []
			reg_inv = []
			for i in range(unrelated_subj_num):
				np.fill_diagonal(fc_corr[i],1)
				inverse.append(np.linalg.inv(fc_corr[i]))
				reg_inv.append(np.linalg.inv(fc_corr[i] + gamma*np.eye(region_num)))
			group_prec = np.mean(inverse, axis=0)
			diff = []
			for i in range(unrelated_subj_num):
				diff.append(np.linalg.norm(reg_inv[i][np.triu_indices(region_num,1)] - group_prec[np.triu_indices(region_num,1)]))
			rmse = np.mean(diff)
			if rmse < min_rmse:
				min_rmse = rmse
				opt_gamma = gamma

		fc_prec = np.zeros([unrelated_subj_num,region_num,region_num])
		for i in range(unrelated_subj_num):
			fc_prec[i] = -np.linalg.inv(fc_corr[i] + opt_gamma*np.eye(region_num))

		mdic = {"C": fc_prec, "gamma": opt_gamma}
		savemat(data_dir + 'FC_' + fc_type + '_unrelated' + str(unrelated_subj_num) + '.mat', mdic)

	else:
		raise ValueError('Please specify fc_type as corr or prec!')

def generate_coupling(unrelated_subj_num, fc_type):

	SC = loadmat(data_dir + 'SC_unrelated' + str(unrelated_subj_num) + '.mat')
	sc = SC['C']

	FC = sio.loadmat(data_dir + 'FC_' + fc_type + '_unrelated' + str(unrelated_subj_num) + '.mat')
	fc = FC['C']

	if sc.shape[0] != fc.shape[0]:
		raise ValueError('The number of subjects in SC doesn\'t match the number of subjects in FC!' )
	elif sc.shape[1] != fc.shape[1]:
		raise ValueError('The number of regions in SC and FC are not the same!')
	elif sc.shape[2] != fc.shape[2]:
		raise ValueError('The number of regions in SC and FC are not the same!')
	elif sc.shape[1] != sc.shape[2]:
		raise ValueError('SC matrices are wrong!')
	elif fc.shape[1] != fc.shape[2]:
		raise ValueError('SC matrices are wrong!')
	else:
		region_num = sc.shape[1]
		regionalcp = np.zeros([unrelated_subj_num, region_num])

		for k in range(unrelated_subj_num):
			for i in range(region_num):

				del_sc = np.delete(sc[k][i], i)
				del_fc = np.delete(fc[k][i], i)

				if (np.count_nonzero(del_fc) == 0) | (np.count_nonzero(del_sc) == 0):
					regionalcp[k][i] = 0
				else:
					regionalcp[k][i], p = stats.pearsonr(del_sc, del_fc)
		 
		mdic = {'scfc_coupling': regionalcp}
		savemat(data_dir + 'regionalcp_' + fc_type + '_unrelated' + str(unrelated_subj_num) + '.mat', mdic)

		mean_regionalcp = np.mean(regionalcp, axis=0)

		np.savetxt(data_dir + 'mean_regionalcp_' + fc_type + '_unrelated' + str(unrelated_subj_num) + '.txt', mean_regionalcp, delimiter=',')


