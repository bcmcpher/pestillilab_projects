#!/bin/bash

echo ''
echo 'Running Left Hemisphere Pairs:'

## for left
set -- `find rois -maxdepth 1 -mindepth 1 -type f -name "lh*" | sed 's#.*/##'`
for a; do
    shift
    for b; do
        printf "%s - %s\n" "$a" "$b"
	touch outs/out_${a}_${b}.txt
    done
done

echo ''
echo 'Running Right Hemisphere Pairs'

## for right
set -- `find rois -maxdepth 1 -mindepth 1 -type f -name "rh*" | sed 's#.*/##'`
for a; do
    shift
    for b; do
        printf "%s - %s\n" "$a" "$b"
	touch outs/out_${a}_${b}.txt
    done
done

echo ''
echo 'Running Left-Right Pairs'

## for intersection
set -- `find rois -maxdepth 1 -mindepth 1 -type f -name "*" | sed 's#.*/##'`
for a; do
    shift
    for b; do
	## if output files of either combo don't exist
	if ([ ! -f outs/*_${a}_${b}.txt ] || [ ! -f outs/*_${b}_${a}.txt ]); then
	    printf "%s - %s\n" "$a" "$b"
	    touch outs/out_${a}_${b}.txt
	fi
    done
done
