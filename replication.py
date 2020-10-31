# the script is for replication of the results in the main text
# sample and out-of-sample reliability
import scipy.io as sio
import numpy as np
import scipy.stats as stats

scfc_dir = './data/scfc_data/'
rep_dir = './data/replication/'
allsubj_dir = './data/allsubjects/'
txt_dir = './txt_file/'
fc_dir = scfc_dir + 'fmriclean_concat_hpf_cc400_gsr_FCcov_997subj/'

rep_subj = np.genfromtxt(rep_dir + 'replicated_subj346.txt', dtype = "int", delimiter = ',', usecols = 0)
all_subj = np.genfromtxt(allsubj_dir + 'subjects_rfMRI_dMRI_complete_997.txt', dtype = "float32", delimiter = ',', usecols = 0)

# load all SC
rep_sc = np.zeros([len(rep_subj),392,392])
sc_all = sio.loadmat(scfc_dir + 'sc_ifod2act_cc400_997subj.mat')
for n in range(len(rep_subj)):
    for m in range(len(all_subj)):
        if all_subj[m] == rep_subj[n]:
            rep_sc[n] = sc_all['sift2volnorm'][0][m]
            break  

# load all FC
rep_fc_corr = np.zeros([len(rep_subj),392,392])
rep_fc_prec = np.zeros([len(rep_subj),392,392])

gamma = 0.3
k = 0
for subj in rep_subj:
	rep_fc_corr[k] = sio.loadmat(fc_dir + str(subj) + '_concat_fmriclean_hpf_cc400_gsr_FCcov.mat')['C']
	rep_fc_prec[k] = np.linalg.inv(rep_fc_corr[k] + gamma*np.eye(392))
	k += 1

# compute scfc couping
rep_cp_corr = np.zeros([len(rep_subj), 392])
rep_cp_prec = np.zeros([len(rep_subj), 392])

for k in range(len(rep_subj)):
    for i in range(392):
        del_sc = np.delete(rep_sc[k][i],i)
        del_fc_corr = np.delete(rep_fc_corr[k][i],i)
        del_fc_prec = np.delete(rep_fc_prec[k][i],i)
        if (np.count_nonzero(del_fc_corr) == 0) | (np.count_nonzero(del_sc) == 0) | (np.count_nonzero(del_fc_prec) == 0):
            rep_cp_corr[k][i] = 0
            rep_cp_prec[k][i] = 0
        else:
            rep_cp_corr[k][i], p = stats.pearsonr(del_sc, del_fc_corr)
            rep_cp_prec[k][i], p = stats.pearsonr(del_sc, -del_fc_prec)

np.savetxt(txt_dir + 'rep_regionalcp_concat_hpf_corr.txt',np.mean(rep_cp_corr, axis=0),delimiter=',')
np.savetxt(txt_dir + 'rep_regionalcp_concat_hpf_prec.txt',np.mean(rep_cp_prec, axis=0),delimiter=',')

# only show regions (p < 0.05)
rp_z = np.arctanh(rep_cp_corr)
meanz_cp = np.mean(rp_z, axis=0)
meanpr_corr = np.tanh(meanz_cp)
n = len(rep_subj)
dist = stats.beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
pvalue_rep_cp_corr = 2*dist.cdf(-abs(meanpr_corr))

meanpr_show_corr = np.zeros(392)
meanpr_show_corr[np.where(pvalue_rep_cp_corr < 0.05)[0]] = meanpr_corr[np.where(pvalue_rep_cp_corr < 0.05)[0]]
np.savetxt(txt_dir + 'rep_regionalcp_sig_concat_hpf_corr.txt',meanpr_show_corr,delimiter=',')

rp_z = np.arctanh(rep_cp_prec)
meanz_cp = np.mean(rp_z, axis=0)
meanpr_prec = np.tanh(meanz_cp)
n = len(rep_subj)
dist = stats.beta(n/2 - 1, n/2 - 1, loc=-1, scale=2)
pvalue_rep_cp_prec = 2*dist.cdf(-abs(meanpr_prec))

meanpr_show_prec = np.zeros(392)
meanpr_show_prec[np.where(pvalue_rep_cp_prec < 0.05)[0]] = meanpr_prec[np.where(pvalue_rep_cp_prec < 0.05)[0]]
np.savetxt(txt_dir + 'rep_regionalcp_sig_concat_hpf_prec.txt',meanpr_show_prec,delimiter=',')


