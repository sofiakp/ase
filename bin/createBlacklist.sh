#!/bin/bash

indiv=$1
vcf=$2

#snps=${MAYAROOT}/rawdata/variants/trio/snps/${indiv}.het.vcf
hg19mask=${MAYAROOT}/rawdata/genomes/masks/wgEncodeHg19ConsensusSignalArtifactRegions.bed
#TODO: personalMask=??
svMask=${MAYAROOT}/rawdata/variants/trio/merged/${indiv}.variants.bed

#cat $hg19mask $svMask | cut -f1-3 | windowBed -a $snps -b stdin -w 100 -v > ${MAYAROOT}/rawdata/variants/trio/snps/${indiv}.het.noblacklist.vcf
cat $hg19mask $svMask | cut -f1-3 | windowBed -a $vcf -b stdin -w 100 -u | cut -f1-2 > ${MAYAROOT}/rawdata/variants/trio/snps/${indiv}.blacklist.txt