#!/bin/bash

## Brent McPherson
## 20160129
## figuring out what went wrong w/ l/r/l2r fiber qsub jobs

## build path
SUBJ=105115
TOPDIR=/N/dc2/projects/lifebid/HCP/Brent
PRJDIR=$TOPDIR/vss-2016
OUTDIR=$PRJDIR/mrtrix
ROIDIR=$OUTDIR/rois
TCKDIR=$OUTDIR/ensemble_wm

set -- `find $ROIDIR -maxdepth 1 -mindepth 1 -type f -name "*.mif" | sed 's#.*/##'`
## set files to be all rois - no path info
for a; do
    shift
    for b; do
        ## a/b are the files w/ extensions
 	## A/B are the files w/o extensions

 	A=`echo $a | sed 's/.mif//g'`
 	B=`echo $b | sed 's/.mif//g'`
 
	for c in SD_STREAM SD_PROB; do
	    for d in 0.25 0.50 1.00 2.00 4.00; do

		CLABEL=`printf "%07g" $COUNT`
		echo tck${CLABEL}_${c}_${d}_${A}_to_${B}.tck >> mrtrix_debug_all.txt

	    done
	done
    done
done
