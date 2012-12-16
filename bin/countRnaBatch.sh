#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Calls qcSample for several samples
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR [Required] Output directory.
   -g FILE [Required] GTF file with gene annotations.
   -e FILE Dir for exon counts 
   -m INT Minimum mapping quality to consider a read. [30]
   -l STR List of samples. If not provided, it will read from STDIN.
   -c     Overwrite 0]
EOF
}

CLEAN=''
INDIR=
OUTDIR=
EXDIR=
GTF=
MQ=30
LIST=
while getopts "hi:o:e:l:g:m:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	e)
	    EXDIR=$OPTARG;;
	g) 
	    GTF=$OPTARG;;
	m)
	    MQ=$OPTARG;;
	l)
	    LIST=$OPTARG;;
	c)
	    CLEAN='-c';;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $OUTDIR || -z $INDIR || -z $GTF ]]; then
    usage; exit 1;
fi
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi
if [[ ! -z $EXDIR && ! -d $EXDIR ]]; then
    mkdir -p $EXDIR
fi
if [ ! -d $INDIR ]; then
    echo 'Indir does not exist' 1>&2; exit 1;
fi
if [ ! -f $GTF ]; then
    echo "GTF does not exits." 1>&2; exit 1;
fi

while read -r sample fq2; do
    if [[ "$sample" =~ ^SNYDER_HG19_ ]]; then
	sample=$sample
    else
	sample=${sample^^} # This will convert to uppercase
	if [[ "$sample" =~ ^[0-9]+ ]]; then # Correct HapMap names
	    sample="GM"$sample
	fi
	sample="SNYDER_HG19_${sample}"
    fi

    infile=${INDIR}/${sample}_reconcile.dedup.bam
    outpref=${OUTDIR}/${sample}
    if [[ ! -z $EXDIR ]]; then
	expref="-e ${EXDIR}/${sample}"
    else
	expref=""
    fi
    if [ ! -f $infile ]; then
	echo "Skipping $sample. Input file is missing." 1>&2; continue;
    fi
    if [[ $fq2 == "NA" ]]; then
	echo "Are you sure you want to count the single-ends?" 1>&2; continue;
	single="-s"
    else
	single=""
    fi
    if [[ $sample =~ RNA || $sample =~ RZ ]]; then
	direction="-d reverse"
    else
	direction="-d no"
    fi
    #if [ ! -s ${outpref}_ambiguous.genecounts ]; then
    bsub -J ${sample}_genes -e /dev/null -o /dev/null -q research-rh6 -W 24:00 -M 8192 -R "rusage[mem=8192]" "${MAYAROOT}/src/bin/countRnaSample.sh -i $infile -o $outpref -g $GTF -m $MQ $expref $single $direction $CLEAN"
    #else
    #	echo "Skipping $sample. Output file exists." 1>&2; continue;
    #fi
done < "${LIST:-/proc/${$}/fd/0}"