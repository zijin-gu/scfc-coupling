# generate surface file and subcortical file
# the input roivalue is the 392 vector saved in a txt file, e.g. roivalue.txt

import numpy as np
import os 
import nibabel as nib
import nibabel.processing

cwd = os.getcwd()
data_dir = cwd + '/data/'

def generate_surface(roivalue):
	#current_path = os.getcwd()
	current_path = cwd
	atlas_dir = current_path + '/atlas_file/'
	image_dir = current_path + '/image_file/'
	txt_dir = current_path + '/herit_results/'

	atlas = nib.load(atlas_dir + 'cc400_roi_atlas.nii.gz').get_fdata()
	nodes = np.unique(atlas)
	nodes = np.delete(nodes,0)
	data = np.zeros(atlas.shape, dtype=np.float32)
	scalar = np.genfromtxt(txt_dir + roivalue + '.txt', dtype = "float32", delimiter = ',', usecols = 0)
	for i,n in enumerate(nodes):
		data[atlas == n] = scalar[i]
	sample_img = nib.load(atlas_dir + 'cc400_roi_atlas.nii.gz')
	save_file = image_dir + roivalue + '_metrics.nii.gz' 
	save_img = nib.Nifti1Image(data, sample_img.affine, sample_img.header)
	save_img.set_data_dtype(data.dtype)
	nib.save(save_img, save_file) 

	filename = save_file
	surf_prefix = filename.strip("cc400_roi_atlas.nii.gz")
	for hemi in ['L', 'R']: 
		os.chdir('/Users/zg243/workbench/bin_macosx64')
		cmd = './wb_command -volume-to-surface-mapping '+  filename + ' /Users/zg243/workbench/HCP_S1200_GroupAvg_v1/S1200.' + hemi + '.midthickness_MSMAll.32k_fs_LR.surf.gii '+ surf_prefix + hemi + '.shape.gii -enclosing'
		os.system(cmd)
		cmd = './wb_command -metric-dilate ' + surf_prefix + hemi + '.shape.gii' + ' /Users/zg243/workbench/HCP_S1200_GroupAvg_v1/S1200.' + hemi + '.midthickness_MSMAll.32k_fs_LR.surf.gii 20 ' + surf_prefix + hemi + '_filled.shape.gii -nearest'
		os.system(cmd)

	os.chdir(current_path)
	os.remove(save_file)
	
def generate_subcortex(roivalue):
	current_path = cwd
	atlas_dir = current_path + '/atlas_file/'
	image_dir = current_path + '/image_file/'
	txt_dir = current_path + '/herit_results/'

	roivol = nib.load(atlas_dir + "cc400_new1mm.nii.gz")
	Vroi = roivol.get_fdata()
	Vnew = np.zeros(Vroi.shape)
	roidata = np.genfromtxt(txt_dir + roivalue + '.txt', dtype = "float32", delimiter = ',', usecols = 0)
	for i,v in enumerate(np.unique(Vroi[Vroi>0])):
		Vnew[Vroi == v] = roidata[i]
	imgnew = nib.Nifti1Image(Vnew, affine = roivol.affine, header = roivol.header)
	nib.save(imgnew, "newdata.nii.gz")

	newdata = nib.load("newdata.nii.gz")
	cc400 = nib.load(atlas_dir + "cc400_new1mm.nii.gz")
	atlas2 = nib.load(atlas_dir + "Atlas_ROIs.2.nii.gz")
	atlas1=nibabel.processing.resample_from_to(atlas2, cc400, order=0)
	V400 = cc400.get_fdata()
	V1 = atlas1.get_fdata()
	subcortvals = np.unique(V400[(V400>0) * (V1>0)])
	V400_subcort = V400 * np.isin(V400,subcortvals)
	cc400_subcort = nib.Nifti1Image(V400_subcort, affine = cc400.affine, header = cc400.header)
	#nib.save(cc400_subcort,"cc400_new1mm_subcort.nii.gz")
	Vnew = newdata.get_fdata()*(V400_subcort>0)
	newdata_subcort=nib.Nifti1Image(Vnew, affine=cc400.affine, header=cc400.header)
	nib.save(newdata_subcort, image_dir + roivalue + "_subcort.nii.gz")
	os.remove("newdata.nii.gz")


def generate_subcortex_fs191(roivalue):
	current_path = cwd
	atlas_dir = current_path + '/atlas_file/'
	image_dir = current_path + '/image_file/'
	txt_dir = current_path + '/glmresults_txt/fs191/'
	roivol = nib.load(atlas_dir + "fs191_sgmfix_dil1_allsubj_mode.nii.gz")
	Vroi = roivol.get_fdata()
	Vroi[Vroi>43]=0
	Vnew = np.zeros(Vroi.shape)
	roidata = np.genfromtxt(txt_dir + roivalue + '.txt', dtype = "float32", delimiter = ',', usecols = 0)
	for i,v in enumerate(np.unique(Vroi[Vroi>0])):
		Vnew[Vroi == v] = roidata[i]
	newdata_subcort = nib.Nifti1Image(Vnew, affine = roivol.affine, header = roivol.header)
	nib.save(newdata_subcort, image_dir + roivalue + "_subcort.nii.gz")

	