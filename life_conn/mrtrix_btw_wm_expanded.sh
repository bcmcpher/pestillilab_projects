#!/bin/bash

## Brent McPherson 
## 20160127
## create hemisphere specific ROIs between 

## module calls
module unload mrtrix/0.3.12
module load mrtrix/0.2.12

## build path
SUBJ=105115
TOPDIR=/N/dc2/projects/lifebid/HCP/Brent
PRJDIR=$TOPDIR/vss-2016
OUTDIR=$PRJDIR/mrtrix
ROIDIR=$OUTDIR/rois
TCKDIR=$OUTDIR/ensemble_tracks2
COUNT=1

set -- `find $ROIDIR -maxdepth 1 -mindepth 1 -type f -name "*.mif" | sed 's#.*/##'`
## set files to be all rois - no path info
for a; do
    shift
    for b; do
        ## a/b are the files w/ extensions
 	## A/B are the files w/o extensions

 	A=`echo $a | sed 's/.mif//g'`
 	B=`echo $b | sed 's/.mif//g'`
 	printf "%s to %s\n" "$A" "$B"
	
	## create specific label for seed in case multiple are running
	SEED=tmp_${A}_${B}.mif
	
	## create seed mask
	mradd $ROIDIR/$a $ROIDIR/$b $ROIDIR/$SEED -quiet

	for c in SD_STREAM SD_PROB; do
	    for d in 0.25 0.50 1.00 2.00 4.00; do

		CLABEL=`printf "%07g" $COUNT`

		echo tck${CLABEL}_${c}_${d}_${A}_to_${B}.tck
		streamtrack $c $OUTDIR/CSD10.mif $TCKDIR/tck${CLABEL}_${c}_${d}_${A}_to_${B}.tck \
                    -seed $ROIDIR/$SEED -mask $OUTDIR/dwi_data_b2000_aligned_trilin_wm_cc.mif \
                    -grad $OUTDIR/dwi_data_b2000_aligned_trilin.b \
                    -include $ROIDIR/$a -include $ROIDIR/$b -number 200 -maxnum 50000 -curvature $d

		let COUNT=COUNT+=1
	    done
	done

	## clear out seed file...
	rm -f $ROIDIR/$SEED

    done
done
