#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Batch convertion from bigWig format to mat format (extractSignal compatible). List of bigWig files is read from stdin.
OPTIONS:
   -h     Show this message and exit
   -o DIR [Required] Output directory.
   -i DIR Input directory. This will be appended to the files read. If this is empty, the lines read are 
          assumed to be paths to files. If this is not empty, the files read are assumed to be basenames.
   -l FILE  File with chromosome lengths [Default: \$GENOMESIZEDIR/hg19.genome].
   -c     Overwrite [0].
EOF
}

INDIR=
OUTDIR=
CHRLEN=${GENOMESIZEDIR}/hg19.genome
CLEAN=0
while getopts "hi:o:l:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	l)
	    CHRLEN=$OPTARG;;
	c)
	    CLEAN=1;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $OUTDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi

JOBGRPID="/bigWig2mat$RANDOM"
bgadd -L 30 $JOBGRPID

while read f; do
    base=$(basename $f)
    if [ -z $INDIR ]; then
	infile=$f
    else
	infile=${INDIR}/$base
    fi
    outpref=${base/.bw/.mat}
    outfile=${OUTDIR}/${outpref}
    if [[ $CLEAN -eq 1 || ! -f $outfile ]]; then
	script=${OUTDIR}/${outpref}.m
	echo "addpath('$MAYAROOT/src/matlab');" > $script
	echo "bigWig2mat('${infile}', '${outfile}', '${CHRLEN}', 1000000, false, '${TMP}');" >> $script
	bsub -g $JOBGRPID -J ${outpref} -o /dev/null -e /dev/null -q research-rh6 -W 24:00 -M 12288 -R "rusage[mem=12288]" "matlab -nodisplay -nosplash < $script"
	#rm $script
    fi
done