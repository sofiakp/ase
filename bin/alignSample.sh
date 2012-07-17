#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Runs BWA on a single sample
OPTIONS:
   -h            Show this message and exit
   --fq1 FILE    [Required]
   --fq2 FILE    [Required] Fastq files
   --bamdir DIR  [Required] Dir where bam files will be written
   --seqpref STR [Required] Prefix for BWA genome index files
   --sample STR  [Required] Sample name. Output will be in <bamdir>/<sample>.bam.
   -c            Overwrite output files [0]
EOF
}

ARGS=`getopt -o "hc" -l "fq1:,fq2:,bamdir:,seqpref:,sample:" -- "$@"`
eval set -- "$ARGS"

CLEAN=0
BAMDIR=
FQ1=
FQ2=
SEQPREF=
SAMPLE=
while [ $# -gt 0 ]; do
    case $1 in
	-h) usage; exit;;
	--fq1) FQ1=$2; shift 2;;
	--fq2) FQ2=$2; shift 2;;
	--bamdir) BAMDIR=$2; shift 2;;
	--seqpref) SEQPREF=$2; shift 2;;
	--sample) SAMPLE=$2; shift 2;;
	-c) CLEAN=1; shift;;
	--) shift; break;;
    esac	    
done

if [ $# -ne 0 ]; then
    usage; exit 1;
fi

if [[ -z $BAMDIR || -z $FQ1 || -z $FQ2 || -z $SEQPREF || -z $SAMPLE ]]; then
    usage; exit 1;
fi

#if [-z $LOGDIR ]
#then
#    LOGDIR=${BAMDIR}/bwaLog
#fi
if [ ! -d $BAMDIR ]; then
    mkdir -p $BAMDIR
fi

# -I not needed if you do SeqPrep first!
# DOUBLE-CHECK THESE PARAMETERS
ALN="-q 20 -t 4"

#### BWA ####
if [[ -s ${FQ1} && -s ${FQ2} ]]; then
    # Check Illumina version
    format=`zcat $FQ1 | python ${MAYAROOT}/src/python/checkFastqFormat.py`
    if [[ $format == '1' ]]; then
	ALN="$ALN -I"
    fi
    
    if [[ $CLEAN -eq 1 || ! -f ${BAMDIR}/${SAMPLE}_1.sai ]]; then
	bwa aln $ALN $SEQPREF ${FQ1} -f ${BAMDIR}/${SAMPLE}_1.sai #2>> ${logfile}
    fi
    if [[ $CLEAN -eq 1 || ! -f ${BAMDIR}/${SAMPLE}_2.sai ]]; then
	bwa aln $ALN $SEQPREF ${FQ2} -f ${BAMDIR}/${SAMPLE}_2.sai #2>> ${logfile}
    fi
else
    echo "Missing or empty FASTQ files. Aborting..." 1>&2 ; exit 1;
fi
    
if [[ -s ${BAMDIR}/${SAMPLE}_1.sai && -s ${BAMDIR}/${SAMPLE}_2.sai ]]; then
    if [[ $CLEAN -eq 1 || ! -f ${BAMDIR}/${SAMPLE}.bam ]]; then
        # Notice the sort by name here...
	head="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:Illumina"
	bwa sampe -r $head $SEQPREF ${BAMDIR}/${SAMPLE}_1.sai ${BAMDIR}/${SAMPLE}_2.sai ${FQ1} ${FQ2} | samtools view -Sbh -t ${SEQPREF}.fa.fai - | samtools sort -n -m 2000000000 - ${BAMDIR}/${SAMPLE} #2>> ${logfile}
	# samtools index ${BAMDIR}/${SAMPLE}.bam
    fi
else
    echo ".sai files missing. Aborting..." 1>&2 ; exit 1;
fi
