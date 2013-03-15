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
   -s DIR [Required] Dir with personal genomes or file to align against if -r is set.
   -l STR List of samples (sample indiv fq1 fq2). If not provided, it will read from STDIN.
   -c     Overwrite [0]
   -r     Align to reference [0]
EOF
}

TOREF=0
CLEAN=''
FQDIR=
BAMDIR=
LIST=
SEQDIR=
while getopts "hf:b:l:crs:" opt
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
	r) 
	    TOREF=1;;
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
if [[ ! -d $SEQDIR && $TOREF -eq 0 ]]; then
    echo 'Seqdir does not exist' 1>&2; exit 1;
elif [[ ! -f $SEQDIR && $TOREF -eq 1 ]]; then
    echo 'Input sequence file does not exist' 1>&2; exit 1;
fi
if [ ! -d $FQDIR ]; then
    echo 'Fqdir does not exist' 1>&2; exit 1;
fi
if [[ $TOREF -eq 1 ]]; then
    SEQDIR=${SEQDIR%.fa}
fi
echo $SEQDIR

while read -r sample indiv fq1 fq2; do
    sample=${sample^^} # This will convert to uppercase
    indiv=${indiv^^}
    if [[ "$indiv" =~ ^[0-9]+$ ]]; then # Correct HapMap names
	indiv="GM"$indiv
	sample="GM"$sample
    fi
    sample="SNYDER_HG19_${sample}"
    if [[ $TOREF -eq 0 ]]; then
	seqpref=${SEQDIR}/${indiv}
	if [[ ! -f ${seqpref}.maternal.fa || ! -f ${seqpref}.paternal.fa ]]; then
	    seqpref=${seqpref}/${indiv}
	    if [[ ! -f ${seqpref}.maternal.fa || ! -f ${seqpref}.paternal.fa ]]; then
		echo "Could not file genome files for $sample. Skipping..." 1>&2 ; continue;
	    fi
	fi
    fi
    fqfile1=${FQDIR}/$(basename $fq1)
    if [[ $fq2 == "NA" ]]; then
	fqfile2='NA'	

    else
	fqfile2=${FQDIR}/$(basename $fq2)
    fi
    if [[ -s $fqfile1 && ( -s $fqfile2 || $fqfile2 == 'NA' ) ]]; then
	if [[ $TOREF -eq 0 ]]; then
	    for par in 'maternal' 'paternal'; do
		if [[ $CLEAN == '-c' || ! -f ${BAMDIR}/${sample}_${par}.bam ]]; then
		    echo $sample
		    bsub -J ${sample}_${par} -o /dev/null -eo ${BAMDIR}/${sample}_${par}_align.err -n 4 -q research-rh6 -W 96:00 -M 12288 -R "rusage[mem=12288]" "${MAYAROOT}/src/bin/alignSample.sh --fq1 $fqfile1 --fq2 $fqfile2 --bamdir $BAMDIR --seqpref ${seqpref}.${par} --sample ${sample}_${par} $CLEAN"
		#else
		#    echo "Skipping ${sample}_${par}. Output file exists" 1>&2
		fi
	    done
	else
	    if [[ $CLEAN == '-c' || ! -f ${BAMDIR}/dedup/${sample}_dedup.bam ]]; then
		bsub -J ${sample} -o /dev/null -eo ${BAMDIR}/${sample}_align.err -n 4 -q research-rh6 -W 96:00 -M 12288 -R "rusage[mem=12288]" "${MAYAROOT}/src/bin/alignSample.sh --fq1 $fqfile1 --fq2 $fqfile2 --bamdir $BAMDIR --seqpref ${SEQDIR} --sample ${sample} $CLEAN -p"
	    else
		echo "Skipping ${sample}. Output file exists" 1>&2
	    fi
	fi
    else
	echo "Could not find fastq files for $sample. Skipping..." 1>&2; continue;
    fi
done < "${LIST:-/proc/${$}/fd/0}"
