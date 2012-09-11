#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Calls alignSample for several samples
OPTIONS:
   -h     Show this message and exit
   -f DIR [Required] Dir with fq files 
   -b DIR [Required] Dir where bam files will be written
   -s DIR [Required] Dir with personal genomes
   -l STR List of samples (sample indiv fq1 fq2). If not provided, it will read from STDIN.
   -c     Overwrite [0]
EOF
}

CLEAN=''
FQDIR=
BAMDIR=
LIST=
SEQDIR=
while getopts "hf:b:l:cs:" opt
do
    case $opt in
	h)
	    usage; exit;;
	f)
	    FQDIR=$OPTARG;;
	b)
	    BAMDIR=$OPTARG;;
	l)
	    LIST=$OPTARG;;
	c) 
	    CLEAN='-c';;
	s) 
	    SEQDIR=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $BAMDIR || -z $FQDIR || -z $SEQDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $BAMDIR ]; then
    mkdir -p $BAMDIR
fi
if [ ! -d $FQDIR ]; then
    echo 'Fqdir does not exist' 1>&2; exit 1;
fi
if [ ! -d $SEQDIR ]; then
    echo 'Seqdir does not exist' 1>&2; exit 1;
fi

while read -r sample indiv fq1 fq2; do
    sample=${sample^^} # This will convert to uppercase
    indiv=${indiv^^}
    if [[ "$indiv" =~ ^[0-9]+$ ]]; then # Correct HapMap names
	indiv="GM"$indiv
	sample="GM"$sample
    fi
    sample="SNYDER_HG19_${sample}"
    seqpref=${SEQDIR}/${indiv}
    if [[ ! -f ${seqpref}.maternal.fa || ! -f ${seqpref}.paternal.fa ]]; then
	seqpref=${SEQDIR}/${indiv}
	if [[ ! -f ${seqpref}.maternal.fa || ! -f ${seqpref}.paternal.fa ]]; then
	    echo "Could not file genome files for $sample. Skipping..." 1>&2 ; continue;
	fi
    fi
    fqfile1=${FQDIR}/$(basename $fq1)
    fqfile2=${FQDIR}/$(basename $fq2)
    if [[ -s $fqfile1 && -s $fqfile2 ]]; then
	for par in 'maternal' 'paternal'; do
	    if [[ $CLEAN == '-c' || ! -f ${BAMDIR}/${sample}_${par}.bam ]]; then
		bsub -J ${sample}_${par} -o /dev/null -eo ${BAMDIR}/${sample}_${par}_align.err -n 4 -q research-rh6 -W 96:00 -M 12288 -R "rusage[mem=12288]" "${MAYAROOT}/src/bin/alignSample.sh --fq1 $fqfile1 --fq2 $fqfile2 --bamdir $BAMDIR --seqpref ${seqpref}.${par} --sample ${sample}_${par} $CLEAN"
	    else
		echo "Skipping ${sample}_${par}. Output file exists" 1>&2
	    fi
	done
    else
	echo "Could not find fastq files for $sample. Skipping..." 1>&2; continue;
    fi
done < "${LIST:-/proc/${$}/fd/0}"
