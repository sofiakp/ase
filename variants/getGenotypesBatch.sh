#!/bin/bash

GDIR=/media/fusion10/work/chromatinVariation/rawdata/genomes/
for f in `ls . | egrep .*snps.vcf`; do
    indiv=${f/.snps.vcf/}
    echo $indiv
    if [ ! -f allNonSan/${indiv}.allNonSan.txt ]; then
	/usr/local/bin/python2.7 ../../supporting/getGenotypes.py $f allNonSan.sites.vcf /media/fusion10/work/chromatinVariation/rawdata/genomes/${indiv}/${indiv}.maternal.fa > allNonSan/${indiv}.allNonSan.txt
    fi
done