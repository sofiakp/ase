#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Counts the number of reads in a set of BAM files.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o FILE [Required] Output directory.
EOF
}

INDIR=
OUTFILE=
while getopts "hi:o:" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || -z $OUTDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $INDIR ]; then
    echo 'Indir does not exist' 1>&2; exit 1;
fi
if [[ ! -d $OUTDIR ]]; then
    mkdir $OUTDIR
fi

#tmpdir=${TMP}/counts_$RANDOM
#for f in `ls $INDIR | egrep SNYDER_HG19_GM18505_POL2_3_reconcile.dedup.maternal.bam$`; do
while read -r f ; do
    if [[ ! ( $f =~ bam$ ) ]]; then
	continue
    fi
    name=`basename $f`
    name=${name/.bam/}
    if [ ! -f ${OUTDIR}/$name ]; then
	bsub -J reads -e /dev/null -o /dev/null "$MAYAROOT/src/bin/simpleCountReads.sh ${INDIR}/${f} $OUTDIR"
    fi
done

#for f in `ls $tmpdir`; do
#    num=`cat ${tmpdir}/$f`
#    echo -e "${f}\t${num}" >> $OUTFILE
#done