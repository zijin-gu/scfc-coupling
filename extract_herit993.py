import scipy.io as sio
import numpy as np
import os
import scipy.stats as stats

cwd = os.getcwd()
data_dir = cwd + '/data/'
herit_dir = data_dir + 'heritability/'
herit_fcdir = herit_dir + 'fmriclean_FCmat_hpf_cc400_gsr_FCcov_997subj_20201022/'

roi_number = 392
all_subj = np.genfromtxt(data_dir + 'subjects_all997.txt', dtype='int')
subjects = np.genfromtxt(herit_dir + 'herit_subjects993.txt', dtype='int')


def extract_herit_sc(sc_type):

    sc = np.zeros([subjects.shape[0], roi_number, roi_number])

    sc_all = sio.loadmat(data_dir + 'sc_ifod2act_cc400_997subj.mat')

    index = all_subj.searchsorted(subjects)
    i = 0
    for idx in index:
        sc[i] = sc_all[sc_type][0][idx]
        i += 1

    mdic = {"C": sc}
    sio.savemat(herit_dir + 'cc400_SC_' + sc_type + '_herit993.mat', mdic)


def extract_herit_fc(gamma):

    fc_corr1 = np.zeros([subjects.shape[0], roi_number, roi_number])
    fc_corr2 = np.zeros([subjects.shape[0], roi_number, roi_number])
    fc_corr3 = np.zeros([subjects.shape[0], roi_number, roi_number])
    fc_corr4 = np.zeros([subjects.shape[0], roi_number, roi_number])

    n = 0
    for subj in subjects:
        fc_corr1[n] = sio.loadmat(herit_fcdir + str(subj) + '_fmriclean/' + str(
            subj) + '_rfMRI_REST1_LR_fmriclean_hpf_cc400_gsr_FCcov.mat')['C']
        fc_corr2[n] = sio.loadmat(herit_fcdir + str(subj) + '_fmriclean/' + str(
            subj) + '_rfMRI_REST1_RL_fmriclean_hpf_cc400_gsr_FCcov.mat')['C']
        fc_corr3[n] = sio.loadmat(herit_fcdir + str(subj) + '_fmriclean/' + str(
            subj) + '_rfMRI_REST2_LR_fmriclean_hpf_cc400_gsr_FCcov.mat')['C']
        fc_corr4[n] = sio.loadmat(herit_fcdir + str(subj) + '_fmriclean/' + str(
            subj) + '_rfMRI_REST2_RL_fmriclean_hpf_cc400_gsr_FCcov.mat')['C']

        n += 1

    fc_prec1 = np.zeros([subjects.shape[0], roi_number, roi_number])
    fc_prec2 = np.zeros([subjects.shape[0], roi_number, roi_number])
    fc_prec3 = np.zeros([subjects.shape[0], roi_number, roi_number])
    fc_prec4 = np.zeros([subjects.shape[0], roi_number, roi_number])

    for i in range(subjects.shape[0]):
        
        np.fill_diagonal(fc_corr1[k], 1)
        np.fill_diagonal(fc_corr2[k], 1)
        np.fill_diagonal(fc_corr3[k], 1)
        np.fill_diagonal(fc_corr4[k], 1)
        fc_prec1[i] = np.linalg.inv(fc_corr1[i] + gamma*np.eye(roi_number))
        fc_prec2[i] = np.linalg.inv(fc_corr2[i] + gamma*np.eye(roi_number))
        fc_prec3[i] = np.linalg.inv(fc_corr3[i] + gamma*np.eye(roi_number))
        fc_prec4[i] = np.linalg.inv(fc_corr4[i] + gamma*np.eye(roi_number))

    mdic = {"C": fc_corr1, "subj": subjects}
    sio.savemat(herit_dir + "cc400_FCcorr_hpf_993subj_session1.mat", mdic)

    mdic = {"C": fc_corr2, "subj": subjects}
    sio.savemat(herit_dir + "cc400_FCcorr_hpf_993subj_session2.mat", mdic)

    mdic = {"C": fc_corr3, "subj": subjects}
    sio.savemat(herit_dir + "cc400_FCcorr_hpf_993subj_session3.mat", mdic)

    mdic = {"C": fc_corr4, "subj": subjects}
    sio.savemat(herit_dir + "cc400_FCcorr_hpf_993subj_session4.mat", mdic)

    mdic = {"C": fc_prec1, "subj": subjects, "gamma": gamma}
    sio.savemat(herit_dir + "cc400_FCprec_hpf_993subj_session1.mat", mdic)

    mdic = {"C": fc_prec2, "subj": subjects, "gamma": gamma}
    sio.savemat(herit_dir + "cc400_FCprec_hpf_993subj_session2.mat", mdic)

    mdic = {"C": fc_prec3, "subj": subjects, "gamma": gamma}
    sio.savemat(herit_dir + "cc400_FCprec_hpf_993subj_session3.mat", mdic)

    mdic = {"C": fc_prec4, "subj": subjects, "gamma": gamma}
    sio.savemat(herit_dir + "cc400_FCprec_hpf_993subj_session4.mat", mdic)


def gen_herit_scfc_couping(sc_type, fc_type, session_id):

    sc = sio.loadmat(herit_dir + 'cc400_SC_' + sc_type + '_herit993.mat')['C']
    fc = sio.loadmat(herit_dir + 'cc400_FC' + fc_type +
                     '_hpf_993subj_session' + str(session_id) + '.mat')['C']
    coupling = np.zeros([subjects.shape[0], roi_number])

    for k in range(subjects.shape[0]):
        for i in range(roi_number):
            del_sc = np.delete(sc[k][i], i)
            del_fc = np.delete(fc[k][i], i)

            if (np.count_nonzero(del_fc) == 0) | (np.count_nonzero(del_sc) == 0):
                coupling[k][i] = 0
            elif fc_type == 'corr':
                coupling[k][i], p = stats.pearsonr(del_sc, del_fc)
            elif fc_type == 'prec':
                coupling[k][i], p = stats.pearsonr(del_sc, -del_fc)

    mdic = {"cp": coupling}
    sio.savemat(herit_dir + 'regionalcp_hpf_' + fc_type +
                '_herit993_session' + str(session_id) + '.mat', mdic)


def standardize(data_type, fc_type=None, session_id=None):

    if data_type != 'cp':

        if data_type == 'sc':
            data = sio.loadmat(herit_dir + 'cc400_SC_sift2volnorm_herit993.mat')['C']

        elif data_type == 'fc':
            data = sio.loadmat(herit_dir + 'cc400_FC' + fc_type +
                               '_hpf_993subj_session' + str(session_id) + '.mat')['C']

        for k in range(data.shape[0]):
            np.fill_diagonal(data[k], 0)

        data_degree = np.sum(abs(data),axis=1)
        data_degree_norm = (data_degree - np.mean(data_degree, axis=1).reshape(
            data.shape[0], 1))/np.std(data_degree, axis=1).reshape(data.shape[0], 1)

        mdic = {"type": data_degree_norm}
        if data_type == 'sc':
            sio.savemat(herit_dir + 'herit2model_' + data_type + '.mat', mdic)
        elif data_type == 'fc':
            sio.savemat(herit_dir + 'herit2model_' + data_type + '_' + fc_type +
                        '_session' + str(session_id) + '.mat', mdic)

    else:
        data = sio.loadmat(herit_dir + 'regionalcp_hpf_' + fc_type +
                           '_herit993_session' + str(session_id) + '.mat')['cp']
        data_norm = (data - np.mean(data, axis=1).reshape(data.shape[0], 1))/np.std(
            data, axis=1).reshape(data.shape[0], 1)

        mdic = {"type": data_norm}
        sio.savemat(herit_dir + 'herit2model_' + data_type + '_' + fc_type +
                    '_session' + str(session_id) + '.mat', mdic)
