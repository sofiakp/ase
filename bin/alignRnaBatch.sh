#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Calls alignSample for several samples
OPTIONS:
   -h     Show this message and exit
   -f DIR [Required] Dir with fq files 
   -b DIR [Required] Dir where tophat files will be written.
   -g FILE [Required] GTF file.
   -m DIR Dir where merged files will be written [Default \$BAMDIR/merged/].
   -s DIR Dir with Bowtie2 indices of personal genomes [Default \$MAYAROOT/rawdata/genomes]
   -t DIR Dir with Bowtie2 indices of personal transcriptomes [Default \$MAYAROOT/rawdata/transcriptomes]
   -l STR List of samples (sample indiv fq1 fq2). If not provided, it will read from STDIN.
   -c     Overwrite [0]
EOF
}

CLEAN=''
FQDIR=
BAMDIR=
MDIR=
GTF=
LIST=
SEQDIR=${MAYAROOT}/rawdata/genomes
TRDIR=${MAYAROOT}/rawdata/transcriptomes
while getopts "hf:b:g:l:cs:t:m:" opt
do
    case $opt in
	h)
	    usage; exit;;
	f)
	    FQDIR=$OPTARG;;
	b)
	    BAMDIR=$OPTARG;;
	g)
	    GTF=$OPTARG;;
	l)
	    LIST=$OPTARG;;
	c) 
	    CLEAN='-c';;
	s) 
	    SEQDIR=$OPTARG;;
	t)
	    TRDIR=$OPTARG;;
	m)
	    MDIR=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $BAMDIR || -z $FQDIR || -z $GTF ]]; then
    usage; exit 1;
fi
if [ ! -d $BAMDIR ]; then
    mkdir -p $BAMDIR
fi
if [ -z $MDIR ]; then
    MDIR=${BAMDIR}/merged
fi
if [ ! -d $MDIR ]; then
    mkdir -p $MDIR
fi
if [ ! -d $FQDIR ]; then
    echo 'Fqdir does not exist' 1>&2; exit 1;
fi
if [ ! -d $SEQDIR ]; then
    echo 'Seqdir does not exist' 1>&2; exit 1;
fi
if [ ! -d $TRDIR ]; then
    mkdir -p $TRDIR
fi
if [ ! -f $GTF ]; then
    echo 'GTF file does not exist' 1>&2; exit 1;
fi

while read -r sample indiv fq1 fq2; do
    sample=${sample^^} # This will convert to uppercase
    indiv=${indiv^^}
    if [[ "$indiv" =~ ^[0-9]+$ ]]; then # Correct HapMap names
	indiv="GM"$indiv
	sample="GM"$sample
    fi
    sample="SNYDER_HG19_${sample}"
    seqpref=${SEQDIR}/${indiv}/${indiv}
    if [[ ! -f ${seqpref}.maternal.1.bt2 || ! -f ${seqpref}.paternal.1.bt2 ]]; then
	echo "Could not file genome files for $sample. Skipping..." 1>&2 ; continue;
    fi
    trpref=${TRDIR}/${indiv}
    if [ ! -d $trpref ]; then
	mkdir $trpref
    fi
    trpref=${trpref}/${indiv}

    fqfile1=${FQDIR}/$(basename $fq1)
    fqfile2=${FQDIR}/$(basename $fq2)
    if [[ -s $fqfile1 && -s $fqfile2 ]]; then
	for par in 'maternal' 'paternal'; do
	    bampref=${BAMDIR}/${sample}_${par}
	    mpref=${MDIR}/${sample}_${par}
	    if [[ $CLEAN == '-c' || ! -f ${mpref}.bam ]]; then
		bsub -J ${sample}_${par} -o /dev/null -eo ${bampref}/align.err -n 4 -q research-rh6 -W 96:00 -M 12288 -R "rusage[mem=12288]" "${MAYAROOT}/src/bin/alignRnaSample.sh --fq1 $fqfile1 --fq2 $fqfile2 --bamdir ${bampref} --mpref $mpref --seqpref ${seqpref}.${par} --trpref ${trpref}.${par} --gtf $GTF $CLEAN"
	    else
		echo "Skipping ${sample}_${par}. Output file exists" 1>&2
	    fi
	done
    else
	echo "Could not find fastq files for $sample. Skipping..." 1>&2; continue;
    fi
done < "${LIST:-/proc/${$}/fd/0}"
