# generate regional txt files
import numpy as np
import scipy.io as sio
import scipy.stats as stats
from scipy.io import savemat

txt_dir = './txt_file/'
unrelated420_dir = './data/unrelated420/'
with open(unrelated420_dir + 'subjects_unrelated420_scfc.txt') as unrelated_subj:
    unrelated_subjects = unrelated_subj.readlines()

unrelated_subjects = [int(x) for x in unrelated_subjects]

# load SC
SC = sio.loadmat(unrelated420_dir + 'cc400_SC_sift2volnorm_unrelated420.mat')
sc = SC['C']

# load FC
# load noconcat FC
FC = []
for fc_type in ['corr', 'prec']:
    FC.append(sio.loadmat(unrelated420_dir + 'cc400_FC' +
                          fc_type + '_nofilt_unrelated420.mat'))
fc_noconcat_nofilt_corr = FC[0]['C']
fc_noconcat_nofilt_prec = FC[1]['C']

# load concat FC
FC_concat = []
for filt_type in ['nofilt', 'bpf', 'hpf']:
    for fc_type in ['corr', 'prec']:
        FC_concat.append(sio.loadmat(unrelated420_dir + 'cc400_FC' +
                                     fc_type + '_concat_' + filt_type + '_unrelated420.mat'))
fc_concat_nofilt_corr = FC_concat[0]['C']
fc_concat_nofilt_prec = FC_concat[1]['C']
fc_concat_bpf_corr = FC_concat[2]['C']
fc_concat_bpf_prec = FC_concat[3]['C']
fc_concat_hpf_corr = FC_concat[4]['C']
fc_concat_hpf_prec = FC_concat[5]['C']

# compute the regional coupling
regionalcp_noconcat_nofilt_corr = np.zeros([len(unrelated_subjects), 392])
regionalcp_noconcat_nofilt_prec = np.zeros([len(unrelated_subjects), 392])
regionalcp_concat_nofilt_corr = np.zeros([len(unrelated_subjects), 392])
regionalcp_concat_nofilt_prec = np.zeros([len(unrelated_subjects), 392])
regionalcp_concat_bpf_corr = np.zeros([len(unrelated_subjects), 392])
regionalcp_concat_bpf_prec = np.zeros([len(unrelated_subjects), 392])
regionalcp_concat_hpf_corr = np.zeros([len(unrelated_subjects), 392])
regionalcp_concat_hpf_prec = np.zeros([len(unrelated_subjects), 392])

for k in range(len(unrelated_subjects)):
    for i in range(392):
        del_sc = np.delete(sc[k][i], i)
        del_fc_noconcat_nofilt_corr = np.delete(
            fc_noconcat_nofilt_corr[k][i], i)
        del_fc_noconcat_nofilt_prec = np.delete(
            fc_noconcat_nofilt_prec[k][i], i)
        del_fc_concat_nofilt_corr = np.delete(fc_concat_nofilt_corr[k][i], i)
        del_fc_concat_nofilt_prec = np.delete(fc_concat_nofilt_prec[k][i], i)
        del_fc_concat_bpf_corr = np.delete(fc_concat_bpf_corr[k][i], i)
        del_fc_concat_bpf_prec = np.delete(fc_concat_bpf_prec[k][i], i)
        del_fc_concat_hpf_corr = np.delete(fc_concat_hpf_corr[k][i], i)
        del_fc_concat_hpf_prec = np.delete(fc_concat_hpf_prec[k][i], i)

        if (np.count_nonzero(del_fc_noconcat_nofilt_corr) == 0) | (np.count_nonzero(del_sc) == 0):
            regionalcp_noconcat_nofilt_corr[k][i] = 0
        else:
            regionalcp_noconcat_nofilt_corr[k][i], p = stats.pearsonr(
                del_sc, del_fc_noconcat_nofilt_corr)

        if (np.count_nonzero(del_fc_noconcat_nofilt_prec) == 0) | (np.count_nonzero(del_sc) == 0):
            regionalcp_noconcat_nofilt_prec[k][i] = 0
        else:
            regionalcp_noconcat_nofilt_prec[k][i], p = stats.pearsonr(
                del_sc, -del_fc_noconcat_nofilt_prec)

        if (np.count_nonzero(del_fc_concat_nofilt_corr) == 0) | (np.count_nonzero(del_sc) == 0):
            regionalcp_concat_nofilt_corr[k][i] = 0
        else:
            regionalcp_concat_nofilt_corr[k][i], p = stats.pearsonr(
                del_sc, del_fc_concat_nofilt_corr)

        if (np.count_nonzero(del_fc_concat_nofilt_prec) == 0) | (np.count_nonzero(del_sc) == 0):
            regionalcp_concat_nofilt_prec[k][i] = 0
        else:
            regionalcp_concat_nofilt_prec[k][i], p = stats.pearsonr(
                del_sc, -del_fc_concat_nofilt_prec)

        if (np.count_nonzero(del_fc_concat_bpf_corr) == 0) | (np.count_nonzero(del_sc) == 0):
            regionalcp_concat_bpf_corr[k][i] = 0
        else:
            regionalcp_concat_bpf_corr[k][i], p = stats.pearsonr(
                del_sc, del_fc_concat_bpf_corr)

        if (np.count_nonzero(del_fc_concat_bpf_prec) == 0) | (np.count_nonzero(del_sc) == 0):
            regionalcp_concat_bpf_prec[k][i] = 0
        else:
            regionalcp_concat_bpf_prec[k][i], p = stats.pearsonr(
                del_sc, -del_fc_concat_bpf_prec)

        if (np.count_nonzero(del_fc_concat_hpf_corr) == 0) | (np.count_nonzero(del_sc) == 0):
            regionalcp_concat_hpf_corr[k][i] = 0
        else:
            regionalcp_concat_hpf_corr[k][i], p = stats.pearsonr(
                del_sc, del_fc_concat_hpf_corr)

        if (np.count_nonzero(del_fc_concat_hpf_prec) == 0) | (np.count_nonzero(del_sc) == 0):
            regionalcp_concat_hpf_prec[k][i] = 0
        else:
            regionalcp_concat_hpf_prec[k][i], p = stats.pearsonr(
                del_sc, -del_fc_concat_hpf_prec)

mdic = {"noconcat_nofilt_corr": regionalcp_noconcat_nofilt_corr,
        "noconcat_nofilt_prec": regionalcp_noconcat_nofilt_prec,
        "concat_nofilt_corr": regionalcp_concat_nofilt_corr,
        "concat_nofilt_prec": regionalcp_concat_nofilt_prec,
        "concat_bpf_corr": regionalcp_concat_bpf_corr,
        "concat_bpf_prec": regionalcp_concat_bpf_prec,
        "concat_hpf_corr": regionalcp_concat_hpf_corr,
        "concat_hpf_prec": regionalcp_concat_hpf_prec}
savemat(unrelated420_dir + "regionalcp_all8types.mat", mdic)

'''
mean_regionalcp_noconcat_nofilt_corr = np.mean(
    regionalcp_noconcat_nofilt_corr, axis=0)
mean_regionalcp_noconcat_nofilt_prec = np.mean(
    regionalcp_noconcat_nofilt_prec, axis=0)
mean_regionalcp_concat_nofilt_corr = np.mean(
    regionalcp_concat_nofilt_corr, axis=0)
mean_regionalcp_concat_nofilt_prec = np.mean(
    regionalcp_concat_nofilt_prec, axis=0)
mean_regionalcp_concat_bpf_corr = np.mean(regionalcp_concat_bpf_corr, axis=0)
mean_regionalcp_concat_bpf_prec = np.mean(regionalcp_concat_bpf_prec, axis=0)
mean_regionalcp_concat_hpf_corr = np.mean(regionalcp_concat_hpf_corr, axis=0)
mean_regionalcp_concat_hpf_prec = np.mean(regionalcp_concat_hpf_prec, axis=0)

np.savetxt(txt_dir + 'regionalcp_noconcat_nofilt_corr.txt',
           mean_regionalcp_noconcat_nofilt_corr, delimiter=',')
np.savetxt(txt_dir + 'regionalcp_noconcat_nofilt_prec.txt',
           mean_regionalcp_noconcat_nofilt_prec, delimiter=',')
np.savetxt(txt_dir + 'regionalcp_concat_nofilt_corr.txt',
           mean_regionalcp_concat_nofilt_corr, delimiter=',')
np.savetxt(txt_dir + 'regionalcp_concat_nofilt_prec.txt',
           mean_regionalcp_concat_nofilt_prec, delimiter=',')
np.savetxt(txt_dir + 'regionalcp_concat_bpf_corr.txt',
           mean_regionalcp_concat_bpf_corr, delimiter=',')
np.savetxt(txt_dir + 'regionalcp_concat_bpf_prec.txt',
           mean_regionalcp_concat_bpf_prec, delimiter=',')
np.savetxt(txt_dir + 'regionalcp_concat_hpf_corr.txt',
           mean_regionalcp_concat_hpf_corr, delimiter=',')
np.savetxt(txt_dir + 'regionalcp_concat_hpf_prec.txt',
           mean_regionalcp_concat_hpf_prec, delimiter=',')
'''
