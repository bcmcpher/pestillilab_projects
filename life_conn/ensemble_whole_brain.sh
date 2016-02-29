#!/bin/bash

module unload mrtrix/0.3.12
module load mrtrix/0.2.12

## Brent McPherson
## 20160217
## run whole brain fiber tracking ensemble parameters
##

## build paths
BASEDIR=/N/dc2/projects/lifebid/HCP/Brent/vss-2016
OUTDIR=$BASEDIR/mrtrix
ROIDIR=$OUTDIR/rois

##
## run all combinations of whole brain tracking
##

streamtrack SD_STREAM $OUTDIR/CSD10.mif $OUTDIR/tck_SD_STREAM_0.25_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 0.25

streamtrack SD_STREAM $OUTDIR/CSD10.mif $OUTDIR/tck_SD_STREAM_0.50_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 0.50

streamtrack SD_STREAM $OUTDIR/CSD10.mif $OUTDIR/tck_SD_STREAM_1.00_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 1.00

streamtrack SD_STREAM $OUTDIR/CSD10.mif $OUTDIR/tck_SD_STREAM_2.00_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 2.00

streamtrack SD_STREAM $OUTDIR/CSD10.mif $OUTDIR/tck_SD_STREAM_4.00_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 4.00

streamtrack SD_PROB $OUTDIR/CSD10.mif $OUTDIR/tck_SD_PROB_0.25_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 0.25

streamtrack SD_PROB $OUTDIR/CSD10.mif $OUTDIR/tck_SD_PROB_0.50_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 0.50

streamtrack SD_PROB $OUTDIR/CSD10.mif $OUTDIR/tck_SD_PROB_1.00_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 1.00

streamtrack SD_PROB $OUTDIR/CSD10.mif $OUTDIR/tck_SD_PROB_2.00_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 2.00

streamtrack SD_PROB $OUTDIR/CSD10.mif $OUTDIR/tck_SD_PROB_4.00_whole-brain.tck \
    -seed $OUTDIR/wm_aseg.mif -mask $OUTDIR/wm_aseg.mif \
    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
    -exclude $OUTDIR/gm_mask_bin.mif -number 20000 -maxnum 2000000 -curvature 4.00


