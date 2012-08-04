#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Reconciles a maternal and a paternal bam into one file.
OPTIONS:
   -h           Show this message and exit
   --indir DIR  [Required] Input directory.
   --outdir DIR [Required] Output directory.
   --sample STR [Required] Sample name. Input files are read from <indir>/<sample>_[mp]aternal.bam and output is written in <outdir>/<sample>_reconcile.bam.
   -c           Overwrite output files [0]
EOF
}

ARGS=`getopt -o "hc" -l "indir:,outdir:,sample:" -- "$@"`
eval set -- "$ARGS"

CLEAN=0
INDIR=
OUTDIR=
SAMPLE=
while [ $# -gt 0 ]; do
    case $1 in
	-h) usage; exit;;
	--indir) INDIR=$2; shift 2;;
	--outdir) OUTDIR=$2; shift 2;;
	--sample) SAMPLE=$2; shift 2;;
	-c) CLEAN=1; shift;;
	--) shift; break;;
	*) usage; exit 1;;
    esac	    
done

if [ $# -ne 0 ]; then
    usage; exit 1;
fi

if [[ -z $INDIR || -z $OUTDIR || -z $SAMPLE ]]; then
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
inpref=${INDIR}/${SAMPLE}
outpref=${OUTDIR}/${SAMPLE}_reconcile
dedpref=${OUTDIR}/dedup/${SAMPLE}_reconcile.dedup

if [[ ( ! -s ${inpref}_maternal.bam ) || ( ! -s ${inpref}_paternal.bam ) ]]; then
    echo "Skipping $SAMPLE. Maternal and/or paternal bam file missing." 1>&2; exit 1;
fi

if [[ $CLEAN -eq 1 || ! -f ${outpref}.bam ]]; then
    ${MAYAROOT}/src/ase_cpp/bin/Ase reconcile rg1=paternal rg2=maternal ${inpref}_paternal.bam ${inpref}_maternal.bam ${tmppref}.bam > ${outpref}.out
    echo "Sorting reconciled file..." 1>&2;
    samtools sort -m 2000000000 ${tmppref}.bam ${outpref}
    echo "Indexing reconciled file..." 1>&2;
    samtools index ${outpref}.bam
fi

if [[ $CLEAN -eq 1 || ! -f ${dedpref}.bam ]]; then
    echo "Removing duplicates..." 1>&2;
    java -Xmx8g -jar ${PICARD}/MarkDuplicates.jar I=${outpref}.bam O=${dedpref}.bam M=${dedpref}.stats AS=true TMP_DIR=${tmpdir} VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 REMOVE_DUPLICATES=true
    samtools index ${dedpref}.bam
fi

if [ -d $tmpdir ]; then
    rm -r $tmpdir
fi
