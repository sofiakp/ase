#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Just strips s/_[12]:N:0:[ACGT]+// from names in a bunch of BAM files.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR Output directory. Must specify unless you set c.
   -c      Just check if files need to be fixed and don't make any changes.
EOF
}

INDIR=
OUTDIR=
CHECK=0
while getopts "hi:o:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	c)
	    CHECK=1;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || ( $CHECK -eq 0 && -z $OUTDIR ) ]]; then
    usage; exit 1;
fi
if [ ! -d $INDIR ]; then
    echo 'Indir does not exist' 1>&2; exit 1;
fi
if [[ $CHECK -eq 0 && ! -d $OUTDIR ]]; then
    mkdir $OUTDIR
fi

for f in `ls $INDIR | egrep .bam$`; do
    infile=${INDIR}/${f}
    outfile=${OUTDIR}/${f}
    num=`samtools view $infile | head -1 | egrep "_[12]:N:0:[ACGT]+" | wc -l`
    if [[ $num -gt 0 ]]; then
	echo "$infile needs fixing." 1>&2
	if [[ $CHECK -eq 0 ]]; then
	    samtools view -h ${infile} |  sed -r 's/_[12]:N:0:[ACGT]+//' | samtools view -Sbh - > ${outfile}
	    if [ -f ${infile}.bai ]; then
		samtools index ${outfile}
	    fi
	fi
    else
	echo "No need to change $infile" 1>&2; continue
    fi    
done
