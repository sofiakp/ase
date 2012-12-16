#!/bin/bash


usage()
{
cat <<EOF
usage: `basename $0` options
Batch convertion from BAM to BED. List of BAM files is read from stdin.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR [Required] Output directory.
   -c     Overwrite [0].
   -p     Merge pairs into fragments [0].
EOF
}

INDIR=
OUTDIR=
BEDPE=0
CLEAN=0
while getopts "hi:o:cp" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	c)
	    CLEAN=1;;
	p)
	    BEDPE=1;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || -z $OUTDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi

JOBGRPID="/bamToBed$RANDOM"
bgadd -L 30 $JOBGRPID

while read f; do
    f=$(basename $f)
    if [[ $BEDPE -eq 0 ]]; then
	outpref=${f/.bam/.bed.gz}
    else
	outpref=${f/.bam/.bedpe.gz}
    fi
    outfile=${OUTDIR}/${outpref}
    echo $outfile
    if [[ $CLEAN -eq 1 || ! -f $outfile ]]; then
	if [[ $BEDPE -eq 0 ]]; then
	    bsub -g $JOBGRPID -J ${outpref} -o /dev/null -e /dev/null "samtools view -bh -q 30 ${INDIR}/$f | bamToBed -i stdin | sort -V | gzip -c > $outfile"
	else
	    # Make sure you select proper pairs or the fragments won't make sense
	    script=${OUTDIR}/${outpref}.sh
	    echo "#!/bin/bash" > $script
	    echo "samtools view -bh -f0x2 ${INDIR}/$f | bamToBed -i stdin -bedpe | awk 'BEGIN{OFS=\"\t\"}{if(\$8 > 30 && \$1 == \$4) print \$1, \$2, \$6, \$7, \$8}' | sort -V | gzip -c > $outfile" >> $script
	    # TODO: If the file is still empty after doing this, then assume it was single end and shift by fragment len.
	    bsub -g $JOBGRPID -J ${outpref} -o /dev/null -e /dev/null < $script
	    rm $script
	fi
    fi
done