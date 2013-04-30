#!/bin/bash

JOBGRPID="/downloadFq$RANDOM"
bgadd -L 10 $JOBGRPID
CLEAN=0
while read -r name file1 file2 ; do
    if [[ ! ( $file1 =~ READ || $file1 == "NA" ) && ( $CLEAN -eq 1 || ! -f `basename $file1` ) ]]; then
	bsub -g $JOBGRPID -J download_${name}_1 -e /dev/null -o /dev/null "wget -nv -N -c ${file1}"
    fi
    if [[ ! ( $file2 =~ READ || $file2 == "NA" ) && ( $CLEAN -eq 1 || ! -f `basename $file2` ) ]]; then
	bsub -g $JOBGRPID -J download_${name}_2 -e /dev/null -o /dev/null "wget -nv -N -c ${file2}"
    fi
done
