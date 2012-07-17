#!/bin/bash

DIR=$1
for file in `ls ${DIR} | egrep "GM192.*[mp]aternal.fa$"`; do
    pref=${file%.fa}
    echo $pref
    bsub -J ${pref} -e ${DIR}/${pref}.idx.err -o /dev/null -q research-rh6 -M 8192 -R "rusage[mem=8192]" "bwa index -a bwtsw -p ${DIR}/${pref} ${DIR}/${file}"
done