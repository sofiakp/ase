#!/bin/bash

INDELS1=${GATK}/bundle/1.5/hg19/1000G_phase1.indels.hg19.vcf
INDELS2=${GATK}/bundle/1.5/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
SNPS=${GATK}/bundle/1.5/hg19/dbsnp_135.hg19.vcf
COV="--standard_covs"

pref=`basename $INFILE`
pref=${pref/.bam}
outpref=${OUTDIR}/${pref}
tmpdir="${TMP}/tmp_${pref}_${RANDOM}_recal"
if [ -d $tmpdir ]; then
    echo "Existing temporary directory! Aborting..." 1>&2; exit 1;
else
    mkdir $tmpdir
fi

if [[ ! -f ${outpref}.realign.intervals ]]; then
    java -Xmx8g -Djava.io.tmpdir=${tmpdir} -jar ${GATK}/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${INFILE} -R ${SEQ} -o ${outpref}.realign.intervals -known ${INDELS1} -known ${INDELS2}
fi
if [[ ! -f ${outpref}.realign.bam ]]; then
    java -Xmx8g -Djava.io.tmpdir=${tmpdir} -jar ${GATK}/GenomeAnalysisTK.jar -T IndelRealigner -I ${INFILE} -R ${SEQ} -targetIntervals ${outpref}.realign.intervals -o ${outpref}.realign.bam -known ${INDELS1} -known ${INDELS2}
fi
if [[ ! -f ${outpref}.realign.cov.cvs ]]; then
    java -Xmx8g -Djava.io.tmpdir=$tmpdir -jar ${GATK}/GenomeAnalysisTK.jar -T CountCovariates -I ${outpref}.realign.bam -R ${SEQ} -knownSites $SNPS -knownSites $INDELS1 -knownSites $INDELS2 $COV -recalFile ${outpref}.realign.cov.cvs -nt 4
fi
if [[ ! -f ${outpref}.recal.bam ]]; then
    java -Xmx8g -Djava.io.tmpdir=$tmpdir -jar ${GATK}/GenomeAnalysisTK.jar -T TableRecalibration -I ${outpref}.realign.bam -R ${SEQ} -o ${outpref}.recal.bam -recalFile ${outpref}.realign.cov.cvs
fi
rm $tmpdir