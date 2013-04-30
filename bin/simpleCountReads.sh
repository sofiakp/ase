#!/bin/bash

# Call through countReads.sh

#INFILE=$1
#OUTDIR=$2

#name=`basename $INFILE`

INDIR=$1
for i in `ls $INDIR | egrep -v "0_recon"| egrep .*dedup.bam$`; do 
    #name=${i/.gz}
    #num=$(zcat ${INDIR}/${i} | wc -l)
    name=${i/.bam/}
    num=$($SAMTOOLS18 view -c ${i}  2>&1)
    echo -e "${name}\t${num}"
done