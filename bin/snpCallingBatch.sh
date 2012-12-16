#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Call variants on a merged set of BAM files.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR [Required] Output directory for merged BAMs.
   -v DIR [Required] Output directory for variant calls. 
   -s STR [Required] Individual.
   -f FILE Fasta file [Default \$MAYAROOT/rawdata/genomes/encodeHg19Male/encodeHg19Male.fa].
   -c     Overwrite [0]
EOF
}

CLEAN=''
INDIR=
OUTDIR=
VARDIR=
INDIV=
FASTA=${MAYAROOT}/rawdata/genomes/encodeHg19Male/encodeHg19Male.fa
while getopts "hi:o:v:s:f:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	v) 
	    VARDIR=$OPTARG;;
	s)
	    INDIV=$OPTARG;;
	f)
	    FASTA=$OPTARG;;
	c) 
	    CLEAN='-c';;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || -z $VARDIR || -z $OUTDIR || -z $INDIV ]]; then
    usage; exit 1;
fi
if [ ! -d $INDIR ]; then
    echo "Input directory does not exist" 1>&2; exit 1;
fi
if [ ! -f $FASTA ]; then
    echo "Input sequence file does not exist" 1>&2; exit 1;
fi
if [[ ! -d $OUTDIR ]]; then
    mkdir -p $OUTDIR
fi
if [[ ! -d $VARDIR ]]; then
    mkdir -p $VARDIR
fi

sample=${INDIV^^} # This will convert to uppercase
if [[ "$sample" =~ ^[0-9]+$ ]]; then # Correct HapMap names
    sample="GM"$sample
fi
sample="SNYDER_HG19_${sample}"

files=
for d in `ls $INDIR | egrep ${sample}.*bam$`; do
    if [[ -z $files ]]; then
	files=${INDIR}/$d
    else
	files="$files ${INDIR}/$d"
    fi
done

if [[ -z $files ]]; then
    echo "No input files." 1>&2; exit 1;
fi
for c in {1..22} X Y; do
    bsub -J ${sample}_snps -o /dev/null -e /dev/null -n 4 -q research-rh6 -W 96:00 -M 12288 -R "rusage[mem=12288]" "${MAYAROOT}/src/bin/snpCalling.sh -i \"$files\" -o $OUTDIR -v $VARDIR -s $sample -r chr${c} -f $FASTA $CLEAN"
done