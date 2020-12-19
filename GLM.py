import scipy.io as sio
import numpy as np
import os
import scipy.stats as stats
from scipy.io import savemat, loadmat
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

cwd = os.getcwd()
# make your own data dir if not exist
if not os.path.exists(cwd + '/data/'):
	os.makedirs(cwd + '/data/')
data_dir = cwd + '/data/'

def generate_covariates(unrelated_subj_num):

	subj_age = np.random.randint(22,37, size=(unrelated_subj_num))

	male_num = np.random.randint(unrelated_subj_num)
	female_num = unrelated_subj_num - male_num
	subj_sex = np.ones(unrelated_subj_num)
	subj_sex[np.random.choice(range(unrelated_subj_num), female_num)] = -1

	subj_cog = np.random.rand(unrelated_subj_num)
	subj_icv = np.random.rand(unrelated_subj_num)
	subj_motion = np.random.rand(unrelated_subj_num)

	np.savetxt(data_dir + 'age_unrelated' + str(unrelated_subj_num) + '.txt', subj_age)
	np.savetxt(data_dir + 'sex_unrelated' + str(unrelated_subj_num) + '.txt', subj_sex)
	np.savetxt(data_dir + 'cog_unrelated' + str(unrelated_subj_num) + '.txt', subj_cog)
	np.savetxt(data_dir + 'icv_unrelated' + str(unrelated_subj_num) + '.txt', subj_icv)
	np.savetxt(data_dir + 'motion_unrelated' + str(unrelated_subj_num) + '.txt', subj_motion)

def glm(unrelated_subj_num, fc_type):
	
	if not os.path.exists(cwd + '/glm_res/'):
		os.makedirs(cwd + '/glm_res/')
	glm_dir = cwd + '/glm_res/'

	# get the age, sex... info
	age = np.genfromtxt(data_dir + 'age_unrelated' + str(unrelated_subj_num) + '.txt')
	sex = np.genfromtxt(data_dir + 'sex_unrelated' + str(unrelated_subj_num) + '.txt')
	cog = np.genfromtxt(data_dir + 'cog_unrelated' + str(unrelated_subj_num) + '.txt')
	ICV = np.genfromtxt(data_dir + 'icv_unrelated' + str(unrelated_subj_num) + '.txt')
	motion = np.genfromtxt(data_dir + 'motion_unrelated' + str(unrelated_subj_num) + '.txt')

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
	
	regionalcp = loadmat(data_dir + 'regionalcp_' + fc_type + '_unrelated' + str(unrelated_subj_num) + '.mat')['scfc_coupling']
	
	roi_number = regionalcp.shape[1]

	regionalcp_z = np.arctanh(regionalcp)

	# GLM
	exog = sm.add_constant(x)
	endog = regionalcp_z 

	p_value = []
	params = []
	conf_int = []
	for i in range(roi_number):
		model = sm.GLM(endog[:,i], exog, sm.families.Gaussian())
		res = model.fit()
		p_value.append(res.pvalues)
		params.append(res.params)
		conf_int.append(res.conf_int())
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
		np.savetxt(glm_dir + 'GLMvar' + str(i) + fc_type + '_allroi.txt', allroi_association, delimiter=',' )

		for j in np.where(p_corrected[0] == True)[0]:
			if params[j,i] > 0:
				sigroi_association[j] = (-np.log(p_corrected[1]))[j]
			else:
				sigroi_association[j] = -(-np.log(p_corrected[1]))[j]
		np.savetxt(glm_dir + 'GLMvar' + str(i) + fc_type + '_sigroi.txt', sigroi_association, delimiter=',' )

		mdic = {'conf_int': conf_int}
		savemat(glm_dir + 'GLM_' + fc_type + '_confint.mat', mdic)
