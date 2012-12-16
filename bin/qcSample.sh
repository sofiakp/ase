#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
SPP qc for a single sample
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input file.
   -o DIR [Required] Output directory.
   -s STR [Required] Stats file.
EOF
}

STATS=
INFILE=
OUTDIR=
while getopts "hi:o:s:" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INFILE=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	s) 
	    STATS=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

pref=`basename $INFILE`
pref=${pref/.bam/}
tmpdir="${TMP}/tmp_${pref}_${RANDOM}_qc"
if [ -d $tmpdir ]; then
    echo "Existing temporary directory! Aborting..." 1>&2; exit 1;
else
    mkdir $tmpdir
fi
tmpfile=${tmpdir}/${pref}.bam 
samtools view -hb -q 20 ${INFILE} > ${tmpfile}
Rscript ${MAYAROOT}/tools/spp_package/run_spp_nodups.R -c=${tmpfile} -savp -odir=$OUTDIR -out=$STATS -tmpdir=$tmpdir
rm -r $tmpdir