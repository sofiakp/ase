#!/bin/bash

sample=$1
bamout=$2
fq=$3
tot=`egrep "Total alignments" $bamout | cut -f2`
pairs=`zcat ${fq} | wc -l | awk '{print $1/4}'`
if [[ tot -eq pairs ]]; then
    echo -e "${sample}\t${tot}\t${pairs}\tOK\n"
else
    echo -e "${sample}\t${tot}\t${pairs}\tWRONG\n"
fi