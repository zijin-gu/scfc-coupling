# extract unrelated420 data from exist all subjects file

import numpy as np
import pandas as pd

alldata_dir = './data/allsubjects/'
unrelated420_dir = './data/unrelated420/'

with open(unrelated420_dir + 'subjects_unrelated420_scfc.txt') as unrelated_subj:
	    unrelated_subjects = unrelated_subj.readlines()
unrelated_subjects = [ int(x) for x in unrelated_subjects ]

def extract_age():
	
	subj_age = pd.read_csv(data_dir + 'brainsizes_subjinfo_voxel_framewisedisp_HCP1200_rfMRI_REST1_RL_with_parcstats_20200114.csv', header=None)
	subj_age = subj_age.values
	all_subj = np.zeros(subj_age.shape[0] - 1)
	for k in range(subj_age.shape[0]-1):
	    all_subj[k] = int(subj_age[k+1][0])
	    
	unrelated_subj = np.zeros(len(unrelated_subjects))
	for k in range(len(unrelated_subjects)):
	    unrelated_subj[k] = int(unrelated_subjects[k])

	age = np.zeros(unrelated_subj.shape[0])

	for i in range(all_subj.shape[0]):
	    for j in range(unrelated_subj.shape[0]):
	        if all_subj[i] == unrelated_subj[j]:
	            age[j] = int(subj_age[i+1][2])       
	
	np.savetxt(unrelated420_dir + 'age_unrelated420.txt', age, delimiter=',')


def extract_ICV():
	
	ICV = open(alldata_dir + 'HCP_997_icv.txt' , "r")
	all_icv = ICV.read().split('\n')
	all_icv.pop(-1)
	size = np.zeros([len(all_icv),2])
	for k in range(len(all_icv)):
	    size[k][0] = int(all_icv[k].split(' ')[0])
	    size[k][1] = float(all_icv[k].split(' ')[1])
	    
	brain_size = np.zeros(len(unrelated_subjects))
	for i in range(len(unrelated_subjects)):
	    for j in range(len(size)):
	        if unrelated_subjects[i] == size[j][0]:
	            brain_size[i] = size[j][1]

	np.savetxt(unrelated420_dir + 'ICV_unrelated420.txt', brain_size, delimiter=',')


def extract_motion():

	motion_dir = alldata_dir + 'HCP_997_mvmtreg_FD/'
	motion = np.zeros(len(unrelated_subjects))
	n = 0
	for subj in unrelated_subjects:
	    fd1lr = open(motion_dir + str(subj) + '_REST1_LR_Framewise_Displacement.txt' , "r")
	    fd1rl = open(motion_dir + str(subj) + '_REST1_RL_Framewise_Displacement.txt' , "r")
	    fd2lr = open(motion_dir + str(subj) + '_REST2_LR_Framewise_Displacement.txt' , "r")
	    fd2rl = open(motion_dir + str(subj) + '_REST2_RL_Framewise_Displacement.txt' , "r")

	    fd1lr = fd1lr.read().split('\n')
	    fd1lr.pop(-1)
	    FD1LR = np.array([float(x) for x in fd1lr])
	    fd1rl = fd1rl.read().split('\n')
	    fd1rl.pop(-1)
	    FD1RL = np.array([float(x) for x in fd1rl])
	    fd2lr = fd2lr.read().split('\n')
	    fd2lr.pop(-1)
	    FD2LR = np.array([float(x) for x in fd2lr])
	    fd2rl = fd2rl.read().split('\n')
	    fd2rl.pop(-1)
	    FD2RL = np.array([float(x) for x in fd2rl])
	    
	    motion[n] = (np.mean(FD1LR) + np.mean(FD1RL) + np.mean(FD2LR) + np.mean(FD2RL))/4
	    n += 1
    
	np.savetxt(unrelated420_dir + 'motion_unrelated420.txt', motion, delimiter=',')


def extract_totalcog():
	subj_cog = pd.read_csv(unrelated420_dir + 'subjects_unrelated420_cognition.csv', header=None)
	subj_cog = subj_cog.values
	Fluid = np.zeros([420,1])
	Early = np.zeros([420,1])
	Total = np.zeros([420,1])
	Crystal = np.zeros([420,1])
	for k in range(420):
	    Fluid[k] = subj_cog[k+1][1]
	    Early[k] = subj_cog[k+1][2]
	    Total[k] = subj_cog[k+1][3]
	    Crystal[k] = subj_cog[k+1][4]
	
	np.savetxt(unrelated420_dir + 'totalcog_unrelated420.txt', Total, delimiter=',')

# men is 1 and women is -1
def extract_sex():
	subj_sex = pd.read_csv(unrelated420_dir + 'sex_420.csv', header=None)
	subj_sex = subj_sex.values
	sex = np.zeros(420)
	for k in range(420):
	    if subj_sex[k+1][1] == 'M':
	        sex[k] = 1
	    else:
	        sex[k] = -1
        
	np.savetxt(unrelated420_dir + 'sex_unrelated420.txt', sex, delimiter=',')

	