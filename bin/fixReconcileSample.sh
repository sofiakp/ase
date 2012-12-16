#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Fixes f**ing naming conventions.
OPTIONS:
   -h           Show this message and exit
   --indir DIR  [Required] Input directory.
   --outdir DIR [Required] Output directory.
   --sample STR [Required] Sample name. Input files are read from <indir>/<sample>_[mp]aternal.bam and output is written in <outdir>/<sample>_reconcile.bam.
   --fq FILE    [Required]
   -c           Overwrite output files [0]
EOF
}

ARGS=`getopt -o "hc" -l "indir:,outdir:,sample:,fq:" -- "$@"`
eval set -- "$ARGS"

CLEAN=0
INDIR=
OUTDIR=
SAMPLE=
FQ=
while [ $# -gt 0 ]; do
    case $1 in
	-h) usage; exit;;
	--indir) INDIR=$2; shift 2;;
	--outdir) OUTDIR=$2; shift 2;;
	--sample) SAMPLE=$2; shift 2;;
	--fq) FQ=$2; shift 2;;
	-c) CLEAN=1; shift;;
	--) shift; break;;
	*) usage; exit 1;;
    esac	    
done

if [ $# -ne 0 ]; then
    usage; exit 1;
fi

if [[ -z $INDIR || -z $OUTDIR || -z $SAMPLE || -z $FQ ]]; then
    usage; exit 1;
fi

if [ ! -d ${OUTDIR}/dedup ]; then
    mkdir -p ${OUTDIR}/dedup
fi

tmpdir="${TMP}/tmp_${SAMPLE}_${RANDOM}"
if [ -d $tmpdir ]; then
    echo "Existing temporary directory! Aborting..." 1>&2; exit 1;
else
    mkdir $tmpdir
fi
tmppref=${tmpdir}/tmp
inpref=${INDIR}/${SAMPLE}_reconcile
outpref=${OUTDIR}/${SAMPLE}_reconcile
dedpref=${OUTDIR}/dedup/${SAMPLE}_reconcile.dedup

if [[ $CLEAN -eq 1 || ! -f ${outpref}.bam ]]; then
    num=`zcat $FQ | head -1 | egrep "_1:N:0:[ACGT]+" | wc -l`
    if [[ $num -gt 0 ]]; then 
	samtools view -h ${inpref}.bam |  sed -r 's/_[12]:N:0:[ACGT]+//' | samtools view -Sbh - > ${outpref}.bam
	samtools index ${outpref}.bam
    else
	cp ${inpref}.bam ${outpref}.bam
	cp ${inpref}.bam.bai ${outpref}.bam.bai
    fi
fi

if [[ $CLEAN -eq 1 || ! -f ${dedpref}.bam ]]; then
    java -Xmx8g -jar ${PICARD}/MarkDuplicates.jar I=${outpref}.bam O=${dedpref}.bam M=${dedpref}.stats AS=true TMP_DIR=${tmpdir} VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 REMOVE_DUPLICATES=true
    samtools index ${dedpref}.bam
fi

if [ -d $tmpdir ]; then
    rm -r $tmpdir
fi
