#!/bin/bash

INDIR=$1
indiv=$2

hg19mask=${MAYAROOT}/rawdata/genomes/masks/wgEncodeHg19ConsensusSignalArtifactRegions.bed
svMask=${INDIR}/merged/${indiv}.variants.bed
snps=${INDIR}/snps/${indiv}.snps.het.vcf # Heterozygous SNPs
clean=${INDIR}/snps/${indiv}.snps.het.noblack.vcf # Het SNPs non overlapping blacklist
#bad=${INDIR}/snps/${indiv}.snps.het.black.txt

if [ ! -d ${INDIR}/masks ]; then
    mkdir ${INDIR}/masks
fi

blackfile=${INDIR}/masks/${indiv}.blacklist.bed
cat $hg19mask $svMask | cut -f1-3 | slopBed -i stdin -b 100 -g /media/fusion10/work/sofiakp/UCSC/hg19/chromLen.txt | mergeBed -i stdin | sort -V > $blackfile
#egrep "#CHROM" $snps > $clean
#intersectBed -a $snps -b $blackfile -v >> $clean
