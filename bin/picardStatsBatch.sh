#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Computes insert size stats for a bunch of samples.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR [Required] Output directory.
   -l STR List of samples. If not provided, it will read from STDIN.
EOF
}

INDIR=
OUTDIR=
LIST=
while getopts "hi:o:l:" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	l)
	    LIST=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $OUTDIR || -z $INDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi
if [ ! -d $INDIR ]; then
    echo 'Indir does not exist' 1>&2; exit 1;
fi

while read -r sample; do
    if [[ "$sample" =~ ^SNYDER_HG19_ ]]; then
	sample=$sample
    else
	sample=${sample^^} # This will convert to uppercase
	if [[ "$sample" =~ ^[0-9]+ ]]; then # Correct HapMap names
	    sample="GM"$sample
	fi
	sample="SNYDER_HG19_${sample}"
    fi
    pref=${sample}_reconcile
    infile=${INDIR}/${pref}.bam
    dupfile=${INDIR}/dedup/${sample}_reconcile.dedup.bam
    if [ ! -f $infile ]; then
	echo "Skipping $sample. Input file is missing." 1>&2; continue;
    fi
    bsub -J ${sample}_stats -e /dev/null -o /dev/null -q research-rh6 -W 24:00 -M 8192 -R "rusage[mem=8192]" "${MAYAROOT}/src/bin/picardStatsSample.sh -i $infile -o $OUTDIR -d $dupfile"
done < "${LIST:-/proc/${$}/fd/0}"