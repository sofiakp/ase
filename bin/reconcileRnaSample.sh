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
   --sample STR [Required] Sample name. Input files are read from <indir>/<sample>_[mp]aternal and output is written in <outdir>/<sample>_reconcile.bam.
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

if [ ! -d ${OUTDIR}/reconcile/dedup ]; then
    mkdir -p ${OUTDIR}/reconcile/dedup
fi

tmpdir="${TMP}/tmp_${SAMPLE}_${RANDOM}"
if [ -d $tmpdir ]; then
    echo "Existing temporary directory! Aborting..." 1>&2; exit 1;
else
    mkdir $tmpdir
fi
tmppref=${tmpdir}/tmp
inpref=${INDIR}/${SAMPLE}
outpref=${OUTDIR}/${SAMPLE}
recpref=${OUTDIR}/reconcile/${SAMPLE}_reconcile
dedpref=${OUTDIR}/reconcile/dedup/${SAMPLE}_reconcile.dedup

if [[ ( ! -s ${inpref}_maternal/accepted_hits.bam ) || ( ! -s ${inpref}_paternal/accepted_hits.bam ) ]]; then
    echo "Skipping $SAMPLE. Maternal and/or paternal bam file missing." 1>&2; exit 1;
fi

if [[ $CLEAN -eq 1 || ! -f ${outpref}_maternal.bam ]]; then
    # Keep only primary pairs
    samtools view -bh -F0x100 -f0x2 ${inpref}_maternal/accepted_hits.bam | samtools sort -n -m 2000000000 - ${outpref}_maternal
fi
if [[ $CLEAN -eq 1 || ! -f ${outpref}_paternal.bam ]]; then
    samtools view -bh -F0x100 -f0x2 ${inpref}_paternal/accepted_hits.bam | samtools sort -n -m 2000000000 - ${outpref}_paternal
fi

#if [[ $CLEAN -eq 1 || ! -f ${recpref}.bam ]]; then
#    ${MAYAROOT}/src/ase_cpp/bin/Ase reconcile rg1=paternal rg2=maternal ${outpref}_paternal.bam ${outpref}_maternal.bam ${tmppref}.bam > ${recpref}.out
#    samtools sort -m 2000000000 ${tmppref}.bam ${recpref}
#    samtools index ${recpref}.bam
#fi

#if [[ $CLEAN -eq 1 || ! -f ${dedpref}.bam ]]; then
#    java -Xmx8g -jar ${PICARD}/MarkDuplicates.jar I=${recpref}.bam O=${dedpref}.bam M=${dedpref}.stats AS=true TMP_DIR=${tmpdir} VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=1000000 REMOVE_DUPLICATES=true
#    samtools index ${dedpref}.bam
#fi

if [ -d $tmpdir ]; then
    rm -r $tmpdir
fi
