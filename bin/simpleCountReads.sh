#!/bin/bash

# Call through countReads.sh

#INFILE=$1
#OUTDIR=$2

#name=`basename $INFILE`
#name=${name/.bam/}
#num=$($SAMTOOLS18 view -c ${INFILE}  2>&1); echo -e "${name}\t${num}" > ${OUTDIR}/$name

INDIR=$1
for i in `ls $INDIR | egrep -v "0_recon"| egrep dedup.bed.gz`; do 
    name=${i/.gz}
    num=$(zcat ${INDIR}/${i} | wc -l)
    echo -e "${name}\t${num}"
done