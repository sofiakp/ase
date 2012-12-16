#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Dedups a bunch of samples.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR [Required] 
   -l STR List of samples. If not provided, it will read from STDIN.
   -c     Overwrite [0]
EOF
}

CLEAN=0
INDIR=
OUTDIR=
LIST=
while getopts "hi:o:l:c" opt
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
	c) 
	    CLEAN=1;;
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

while read -r sample ; do
    if [[ "$sample" =~ ^SNYDER_HG19_ ]]; then
	sample=$sample
    else
	sample=${sample^^} # This will convert to uppercase
	if [[ "$sample" =~ ^[0-9]+ ]]; then # Correct HapMap names
	    sample="GM"$sample
	fi
	sample="SNYDER_HG19_${sample}"
    fi
    inpref=${INDIR}/${sample}
    outpref=${OUTDIR}/${sample}_dedup

    if [[ ! -s ${inpref}.bam ]]; then
	echo "Skipping $sample. Input file missing." 1>&2; continue;
    fi
    if [[ $CLEAN -eq 1 || ! -f ${outpref}.bam ]]; then
	bsub -J ${sample}_dedup -e /dev/null -o /dev/null -n 1 -q research-rh6 -W 24:00 -M 16384 -R "rusage[mem=16384]" "java -Xmx8g -jar ${PICARD}/MarkDuplicates.jar I=${inpref}.bam O=${outpref}.bam M=${outpref}.stats AS=true TMP_DIR=${OUTDIR} VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 REMOVE_DUPLICATES=true; samtools index ${outpref}.bam"
    else
	echo "Skipping $sample. Output file exists." 1>&2; continue;
    fi
done < "${LIST:-/proc/${$}/fd/0}"
