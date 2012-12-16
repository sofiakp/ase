#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Picard stats for a bam file.
OPTIONS:
   -h      Show this message and exit
   -i FILE [Required] Input file.
   -d FILE [Required] Input file with duplicates removed.
   -o DIR  [Required] Output directory.
EOF
}

INFILE=
DUPFILE=
OUTDIR=
while getopts "hi:o:d:" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INFILE=$OPTARG;;
	d)
	    DUPFILE=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

pref=`basename $INFILE`
pref=${pref/.bam/}
outpref=${OUTDIR}/${pref}

#if [[ ! -f ${outpref}_qual.txt ]]; then
#    java -Xmx8g -jar ${PICARD}/QualityScoreDistribution.jar I=$INFILE O=${outpref}_qual.txt CHART=${outpref}_qual.pdf ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
#fi 
if [[ ! -f ${outpref}_insSize.txt ]]; then
    java -Xmx8g -jar ${PICARD}/CollectInsertSizeMetrics.jar I=$DUPFILE O=${outpref}_insSize.txt H=${outpref}_insSize.pdf ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
fi 
