import scipy.io as sio
import numpy as np
import os
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

cwd = os.getcwd()
data_dir = cwd + '/data/'
covariate_dir = data_dir + 'covariates/'

if not os.path.exists(cwd + 'glmresults_txt/'):
    os.makedirs(cwd + 'glmresults_txt/')

glm_dir = cwd + '/glmresults_txt/'


roi_number = 392

def glm(fc_type, filter_type):
	
	# get the age, sex... info
	subj_age = np.genfromtxt(covariate_dir + 'age_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)
	subj_sex = np.genfromtxt(covariate_dir + 'sex_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)
	subj_cog = np.genfromtxt(covariate_dir + 'totalcog_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)
	subj_ICV = np.genfromtxt(covariate_dir + 'ICV_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)
	subj_motion = np.genfromtxt(covariate_dir + 'motion_unrelated420.txt', dtype = "float32", delimiter = ',', usecols = 0)

	# remove the subjects who don't have cognitive data
	nan_idx = list(np.argwhere(np.isnan(subj_cog)).reshape(np.argwhere(np.isnan(subj_cog)).shape[0]))

	age = np.delete(subj_age, nan_idx)
	sex = np.delete(subj_sex, nan_idx)
	cog = np.delete(subj_cog, nan_idx)
	motion = np.delete(subj_motion, nan_idx)
	ICV = np.delete(subj_ICV, nan_idx)

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
	x = np.concatenate((cen_age.reshape(-1,1), sex.reshape(-1,1), cen_cog.reshape(-1,1), cen_ICV.reshape(-1,1), cen_motion.reshape(-1,1), agecog.reshape(-1,1),sexcog.reshape(-1,1), ICVmotion.reshape(-1,1)), axis=1)
	
	subj_regionalcp = sio.loadmat(data_dir + 'regionalcp_' + filter_type + '_' + fc_type + '.mat')['scfc_coupling']
	regionalcp = np.delete(subj_regionalcp, nan_idx, axis = 0)
	regionalcp_z = np.arctanh(regionalcp)

	# GLM
	exog = sm.add_constant(x)
	endog = regionalcp_z 

	p_value = []
	params = []
	for i in range(roi_number):
		model = sm.GLM(endog[:,i], exog, sm.families.Gaussian())
		res = model.fit()
		p_value.append(res.pvalues)
		params.append(res.params)
	p_value = np.array(p_value)
	params = np.array(params)

	for i in range(1,exog.shape[1]):
		p_corrected = multipletests(p_value[:,i], alpha = 0.05, method = 'fdr_bh')

		allroi_association = np.zeros(roi_number)
		sigroi_association = np.zeros(roi_number)
		for j in range(roi_number):
		    if params[j,i] > 0:
		        allroi_association[j] = (-np.log(p_corrected[1]))[j]
		    else:
		        allroi_association[j] = -(-np.log(p_corrected[1]))[j]
		np.savetxt(glm_dir + 'GLMvar' + str(i) + '_' + filter_type + '_' + fc_type + '_allroi.txt', allroi_association, delimiter=',' )

		for j in np.where(p_corrected[0] == True)[0]:
		    if params[j,i] > 0:
		        sigroi_association[j] = (-np.log(p_corrected[1]))[j]
		    else:
		        sigroi_association[j] = -(-np.log(p_corrected[1]))[j]
		np.savetxt(glm_dir + 'GLMvar' + str(i) + '_' + filter_type + '_' + fc_type + '_sigroi.txt', sigroi_association, delimiter=',' )
