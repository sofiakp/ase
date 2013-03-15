#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Merges BED files. Reads(from STDIN) a metafile with replicates in the same format used for signal track generation. The first column of the file has replicates separated by commas. These are merged.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR [Required] Output directory.
   -c     Overwrite.
EOF
}

INDIR=
OUTDIR=
CLEAN=0
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
	    CLEAN=1;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || -z $OUTDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $INDIR ]; then
    echo "Input directory does not exist" 1>&2; exit 1;
fi
if [[ ! -d $OUTDIR ]]; then
    mkdir -p $OUTDIR
fi

JOBGRPID="/mergeChromatinSamples$RANDOM"
bgadd -L 30 $JOBGRPID 

while read inputline; do
    chip=$(echo ${inputline} | awk '{print $1}' | sed -r 's/\.bam/.bed.gz/g') # extract first column as ChIP file name, replace .bam with .bed.gz
    nchip=$(echo ${chip} | sed -r 's/;/\n/g' | wc -l)
    chipstub=$(echo ${inputline} | awk '{print $3}')
    outchip=${INDIR}/${chipstub}_0_dedup.bed.gz #0_reconcile.dedup.bed.gz
    if [[ $nchip -gt 1  && ( ! -f $outchip || $CLEAN -eq 1 ) ]]; then
	#logfile=${INDIR}/${chipstub}_0_reconcile.dedup.files.txt # Write which files were used in the merged file
	#echo -e "${chip}\t${chipstub}_0_reconcile.dedup.bed.gz" > $logfile
	chip=$(echo ${chip} | sed -r 's/;/\n/g' | xargs -I fname find "${INDIR}" -name fname -printf "%p ") # separate file names with space
	bsub -q research-rh6 -g "${JOBGRPID}" -J "${chipstub}" -W 24:00 -o /dev/null -e /dev/null "zcat $chip | sort -V | gzip -c > $outchip"
    else
	echo "No need to merge $chip" &>2
    fi
done
