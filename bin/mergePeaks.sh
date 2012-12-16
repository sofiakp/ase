#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Merges narrowPeak files for the given marks.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR [Required] Output directory.
   -m FILE Metafile
EOF
}

INDIR=
OUTDIR=
while getopts "hi:o:m:" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	m)
	    META=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || -z $OUTDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $INDIR ]; then echo "Input directory does not exist" 1>&2; exit 1; fi
if [ ! -f $META ]; then echo "Metafile does not exist" 1>&2; exit 1; fi
if [[ ! -d $OUTDIR ]]; then
    mkdir -p $OUTDIR
fi

while read mark; do
    files=$(egrep $mark $META | awk '{print $1"_VS_"$2"_peaks.encodePeak*"}' | xargs -I fname find $INDIR -name fname -printf "%p ")
    logfile="${OUTDIR}/SNYDER_HG19_${mark}_merged.encodePeak.files.txt"
    outfile="${OUTDIR}/SNYDER_HG19_${mark}_merged.bed.gz"
    echo $files > $logfile
    bsub -q research-rh6 -J ${mark}_merge -W 24:00 -o /dev/null -e /dev/null "zcat $files | cut -f1-3 | sort -V | mergeBed -i stdin | gzip -c > $outfile"
done
