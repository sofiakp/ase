#!/bin/bash

# Call through countReads.sh

INFILE=$1
OUTDIR=$2

name=`basename $INFILE`
name=${name/.bam/}
num=$($SAMTOOLS18 view -c ${INFILE}  2>&1); echo -e "${name}\t${num}" > ${OUTDIR}/$name