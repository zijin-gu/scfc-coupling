# GLM analysis of regional coupling

import numpy as np
import scipy.io as sio
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

unrelated420_dir = './data/unrelated420/'
txt_dir = './txt_file/'
# get the age, sex... info
subj_age = np.genfromtxt(unrelated420_dir + 'age_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)
subj_sex = np.genfromtxt(unrelated420_dir + 'sex_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)
subj_cog = np.genfromtxt(unrelated420_dir + 'totalcog_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)
subj_ICV = np.genfromtxt(unrelated420_dir + 'ICV_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)
subj_motion = np.genfromtxt(unrelated420_dir + 'motion_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)

# remove the subjects who don't have cognitive data
age = np.delete(subj_age, [19, 256,258,317,391])
sex = np.delete(subj_sex, [19, 256,258,317,391])
cog = np.delete(subj_cog, [19, 256,258,317,391])
motion = np.delete(subj_motion, [19, 256,258,317,391])
ICV = np.delete(subj_ICV, [19, 256,258,317,391])

# center all variables
cen_age = age - np.mean(age)
cen_cog = cog - np.mean(cog)
cen_ICV = ICV - np.mean(ICV)
cen_motion = motion - np.mean(motion)

# compute the interactions
sexcog = sex * cen_cog
agecog = cen_age * cen_cog
ICVmotion = cen_ICV * cen_motion

# concatenate all the covariates
x = np.concatenate((cen_age.reshape(-1,1), sex.reshape(-1,1), cen_cog.reshape(-1,1), agecog.reshape(-1,1),sexcog.reshape(-1,1),
                   cen_ICV.reshape(-1,1), cen_motion.reshape(-1,1), ICVmotion.reshape(-1,1)), axis=1)

# load the regional coupling data
regional_cp = sio.loadmat(unrelated420_dir + 'regionalcp_all8types.mat')
regionalcp_noconcat_nofilt_corr = np.delete(regional_cp['noconcat_nofilt_corr'], [19, 256,258,317,391], axis=0)
regionalcp_noconcat_nofilt_prec = np.delete(regional_cp['noconcat_nofilt_prec'], [19, 256,258,317,391], axis=0)
regionalcp_concat_nofilt_corr = np.delete(regional_cp['concat_nofilt_corr'], [19, 256,258,317,391], axis=0)
regionalcp_concat_nofilt_prec = np.delete(regional_cp['concat_nofilt_prec'], [19, 256,258,317,391], axis=0)
regionalcp_concat_bpf_corr = np.delete(regional_cp['concat_bpf_corr'], [19, 256,258,317,391], axis=0)
regionalcp_concat_bpf_prec = np.delete(regional_cp['concat_bpf_prec'], [19, 256,258,317,391], axis=0)
regionalcp_concat_hpf_corr = np.delete(regional_cp['concat_hpf_corr'], [19, 256,258,317,391], axis=0)
regionalcp_concat_hpf_prec = np.delete(regional_cp['concat_hpf_prec'], [19, 256,258,317,391], axis=0)


# z-transform
regionalcp_noconcat_nofilt_corr_z = np.arctanh(regionalcp_noconcat_nofilt_corr)
regionalcp_noconcat_nofilt_prec_z = np.arctanh(regionalcp_noconcat_nofilt_prec)
regionalcp_concat_nofilt_corr_z = np.arctanh(regionalcp_concat_nofilt_corr)
regionalcp_concat_nofilt_prec_z = np.arctanh(regionalcp_concat_nofilt_prec)
regionalcp_concat_bpf_corr_z = np.arctanh(regionalcp_concat_bpf_corr)
regionalcp_concat_bpf_prec_z = np.arctanh(regionalcp_concat_bpf_prec)
regionalcp_concat_hpf_corr_z = np.arctanh(regionalcp_concat_hpf_corr)
regionalcp_concat_hpf_prec_z = np.arctanh(regionalcp_concat_hpf_prec)


# GLM
exog = sm.add_constant(x)
# noconcat
for fc_type in ['corr','prec']:
	endog = eval('regionalcp_noconcat_nofilt_' + fc_type + '_z')

	p_value = []
	params = []
	for i in range(392):
		model = sm.GLM(endog[:,i], exog, sm.families.Gaussian())
		res = model.fit()
		p_value.append(res.pvalues)
		params.append(res.params)
	p_value = np.array(p_value)
	params = np.array(params)

	for i in range(1,4):
		p_corrected = multipletests(p_value[:,i], alpha=0.05, method='fdr_bh')

		allroi_association = np.zeros(392)
		sigroi_association = np.zeros(392)
		for j in range(392):
		    if params[j,i] > 0:
		        allroi_association[j] = (-np.log(p_corrected[1]))[j]
		    else:
		        allroi_association[j] = -(-np.log(p_corrected[1]))[j]
		np.savetxt(txt_dir + 'GLMvar' + str(i) + '_noconcat_nofilt_' + fc_type + '_allroi.txt', allroi_association, delimiter=',' )

		for j in np.where(p_corrected[0] == True)[0]:
		    if params[j,i] > 0:
		        sigroi_association[j] = (-np.log(p_corrected[1]))[j]
		    else:
		        sigroi_association[j] = -(-np.log(p_corrected[1]))[j]
		np.savetxt(txt_dir + 'GLMvar' + str(i) + '_noconcat_nofilt_' + fc_type + '_sigroi.txt', sigroi_association, delimiter=',' )

# concat
for filter_type in ['nofilt','bpf','hpf']:
	for fc_type in ['corr','prec']:
		endog = eval('regionalcp_concat_' + filter_type + '_' + fc_type + '_z')

		p_value = []
		params = []
		for i in range(392):
			model = sm.GLM(endog[:,i], exog, sm.families.Gaussian())
			res = model.fit()
			p_value.append(res.pvalues)
			params.append(res.params)
		p_value = np.array(p_value)
		params = np.array(params)

		for i in range(1,4):
			p_corrected = multipletests(p_value[:,i], alpha=0.05, method='fdr_bh')

			allroi_association = np.zeros(392)
			sigroi_association = np.zeros(392)
			for j in range(392):
			    if params[j,i] > 0:
			        allroi_association[j] = (-np.log(p_corrected[1]))[j]
			    else:
			        allroi_association[j] = -(-np.log(p_corrected[1]))[j]
			np.savetxt(txt_dir + 'GLMvar' + str(i) + '_concat_' + filter_type + '_' + fc_type + '_allroi.txt', allroi_association, delimiter=',' )

			for j in np.where(p_corrected[0] == True)[0]:
			    if params[j,i] > 0:
			        sigroi_association[j] = (-np.log(p_corrected[1]))[j]
			    else:
			        sigroi_association[j] = -(-np.log(p_corrected[1]))[j]
			np.savetxt(txt_dir + 'GLMvar' + str(i) + '_concat_' + filter_type + '_' + fc_type + '_sigroi.txt', sigroi_association, delimiter=',' )

