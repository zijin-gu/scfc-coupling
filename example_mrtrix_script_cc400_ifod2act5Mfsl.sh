SUBJECT=100206

INPUTDIR=HCPDATADIR/${SUBJECT}
diffdir=${INPUTDIR}/T1w/Diffusion
regimg=${diffdir}/nodif_brain_mask.nii.gz

#Convert data and compute CSD
mrconvert $diffdir/data.nii.gz DWI.mif -fslgrad $diffdir/bvecs $diffdir/bvals -datatype float32 -stride 0,0,0,1 -quiet -force
dwibiascorrect DWI.mif DWI_bc.mif -ants -bias DWI_biasfield.mif -mask $diffdir/nodif_brain_mask.nii.gz  -force
dwi2response dhollander DWI_bc.mif RF_wm_dhollander.txt RF_gm_dhollander.txt RF_csf_dhollander.txt -voxels RF_voxels_dhollander.mif  -force

#Make 5TT tissue model using FSL FAST, resampled to diffusion volume space
5ttgen fsl ${INPUTDIR}/T1w/T1w_acpc_dc_restore.nii.gz 5TT_fsl_T1space.mif -mask ${INPUTDIR}/T1w/T1w_acpc_dc_restore_brain.nii.gz -nocrop -force
mrconvert 5TT_fsl_T1space.mif 5TT_fsl_T1space.nii.gz -force
rm -f 5TT_fsl_T1space.mif
$FSLDIR/bin/applywarp --rel --interp=trilinear -i 5TT_fsl_T1space.nii.gz -r ${regimg} -o 5TT_fsl_dwi.nii.gz

#Warp cc400 atlas from MNI space to diffusion volume space
atlasfile="cc400_new1mm_seq392.nii.gz"
roiname="cc400"
roifile="cc400_dwi.nii.gz"
$FSLDIR/bin/applywarp -i ${atlasfile} -r ${regimg} -w ${INPUTDIR}/MNINonLinear/xfms/standard2acpc_dc.nii.gz -o ${roifile} --interp=nn

################################
algo=iFOD2
algo_str="ifod2act5Mfsl"

tckfile=${SUBJECT}_CSD_${algo_str}_cut05_5M.tck
siftfile=${tckfile/.tck/_sift2.txt}
outcsv=${tckfile/.tck/""}_${roiname}_connectome_sift2_volnorm.csv

#Run tractography with 5TT tissue model, seeding dynamically until 5M tracts
tckgen RF_wm_dhollander.mif ${tckfile} -algorithm ${algo} -seed_dynamic RF_wm_dhollander.mif -select 5M -maxlength 300 -cutoff 0.05 -act 5TT_fsl_dwi.nii.gz -force

#Compute SIFT2 tract weights
tcksift2 ${tckfile} RF_wm_dhollander.mif ${siftfile} -fd_thresh 0.05 -force

#note: default endpoint assignment is to search within 4mm of unassigned endpoints for the nearest ROI (-assignment_radial_search 4)
tck2connectome -symmetric -scale_invnodevol ${tckfile} ${roifile} ${outcsv} -tck_weights_in ${siftfile} -force
