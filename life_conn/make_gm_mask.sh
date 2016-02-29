#!/bin/bash

module unload mrtrix/0.3.12
module load mrtrix/0.2.12
module load fsl

## Brent McPherson
## 20160217
## combine region ROIs into gray matter mask for whole brain tracking exclusion mask
##

$BASEDIR=/N/dc2/projects/lifebid/HCP/Brent/vss-2016
$MRTDIR=$BASEDIR/mrtrix

$ROIDIR=$MRXDIR/rois
$OUTDIR=$MRTDIR

## combine all the ROIs into a single mask
mradd $ROIDIR/*.mif $OUTDIR/gm_mask.mif

## convert mask to a nifti
mrconvert $OUTDIR/gm_mask.mif gm_mask.nii.gz

## convert mask to a binary image in nifti format
fslmaths gm_mask.nii.gz -bin gm_mask_bin.nii.gz

## convert binary mask back to .mif format
mrconvert gm_mask_bin.nii.gz gm_mask_bin.mif

