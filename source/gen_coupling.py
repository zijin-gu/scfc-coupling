import scipy.io as sio
import numpy as np
import os
import scipy.stats as stats

cwd = os.getcwd()
data_dir = cwd + '/data/'

if not os.path.exists(cwd + '/regionalcp_txt/'):
    os.makedirs(cwd + 'regionalcp_txt/')
txt_dir = cwd + 'regionalcp_txt/'

all_subj = np.genfromtxt(data_dir + 'subjects_all997.txt', dtype='int')
unrelated_subj = np.genfromtxt(
    data_dir + 'subjects_unrelated420.txt', dtype='int')

roi_number = 392


def gen_scfc_coupling(sc_type, fc_type, filter_type):

    # load SC
    SC = sio.loadmat(data_dir + 'cc400_SC_' + sc_type + '_unrelated420.mat')
    sc = SC['C']

    FC = sio.loadmat(data_dir + 'cc400_FC' + fc_type +
                     '_concat_' + filter_type + '_unrelated420.mat')
    fc = FC['C']

    regionalcp = np.zeros([unrelated_subj.shape[0], roi_number])

    for k in range(unrelated_subj.shape[0]):
        for i in range(roi_number):

            del_sc = np.delete(sc[k][i], i)
            del_fc = np.delete(fc[k][i], i)

            if (np.count_nonzero(del_fc) == 0) | (np.count_nonzero(del_sc) == 0):
            	regionalcp[k][i] = 0
            elif fc_type == 'corr':
                regionalcp[k][i], p = stats.pearsonr(del_sc, del_fc)
            elif fc_type == 'prec':
                regionalcp[k][i], p = stats.pearsonr(del_sc, -del_fc)

    mdic = {'scfc_coupling': regionalcp}
    sio.savemat(data_dir + 'regionalcp_' + filter_type +
                '_' + fc_type + '.mat', mdic)

    mean_regionalcp = np.mean(regionalcp, axis=0)

    np.savetxt(txt_dir + 'regionalcp_' + filter_type + '_' +
               fc_type + '.txt', mean_regionalcp, delimiter=',')
