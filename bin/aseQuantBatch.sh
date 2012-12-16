#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Counts reads with the REF and ALT alleles at the positions specified.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR [Required] Output directory
   -v VCF Vcf file with positions where allele counts will be computed [Default: \$MAYAROOT/rawdata/variants/all/snps/<indiv>.snps.vcf]
   -k FILE Mask file (R data file) with information about the variants [default: same as VCF with suffix r].
   -l STR List of input samples. If not provided, it will read from STDIN.
   -m INT Minimum mapping quality to consider a read. [30]
   -b INT Minimum base quality to consider a position in a read. [20]
   -c     Overwrite [0]
EOF
}

CLEAN=''
INDIR=
OUTDIR=
MASK=
LIST=
VCF=
MQ=30
BQ=20
while getopts "hi:o:l:v:m:b:k:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	l)
	    LIST=$OPTARG;;
	v) 
	    VCF=$OPTARG;;
	k) 
	    MASK=$OPTARG;;
	m)
	    MQ=$OPTARG;;
	b)
	    BQ=$OPTARG;;
	c) 
	    CLEAN='-c';;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $OUTDIR || -z $INDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi
if [ ! -d $INDIR ]; then
    echo 'Input directory does not exist' 1>&2; exit 1;
fi

while read -r sample ; do
    #pref=`basename $infile`
    #pref=${pref/.bam/}

    if [[ "$sample" =~ ^SNYDER_HG19_ ]]; then
	slist=(`echo $sample | tr "_" "\n"`)
	indiv=${slist[2]}
    else
	sample=${sample^^} # This will convert to uppercase
	if [[ "$sample" =~ ^[0-9]+ ]]; then # Correct HapMap names
	    sample="GM"$sample
	fi
	slist=(`echo $sample | tr "_" "\n"`)
	indiv=${slist[0]}
	echo $indiv
	sample="SNYDER_HG19_${sample}"
    fi
    infile=${INDIR}/${sample}_reconcile.dedup.bam
    outpref=${OUTDIR}/${sample}
    if [[ -z $VCF ]]; then
	ivcf=${MAYAROOT}/rawdata/variants/all/snps/allNonSan/${indiv}.allNonSan.txt #snps.vcf
    else
	ivcf=$VCF
    fi
    if [ ! -f $ivcf ]; then
	echo 'VCF file does not exist' 1>&2; continue;
    fi
    if [[ -z $MASK ]]; then
	imask=${ivcf/.vcf/.r}
    else
	imask=$MASK
    fi
    #if [[ ! -f $imask  ]]; then
    #	echo "Mask file does not exist." 1>&2; continue;
    #fi
	
    if [[ ! -s ${infile} || ! -s ${infile}.bai ]]; then
	echo "Skipping $sample. Input bam file missing." 1>&2; continue;
    fi

    if [[ $CLEAN == "-c" || ! -f ${outpref}.counts ]]; then
	bsub -J ${sample}_count -e /dev/null -o /dev/null -q research-rh6 -W 24:00 -M 8192 -R "rusage[mem=8192]" "${MAYAROOT}/src/ase_cpp/bin/Ase asequantmultirg -m $MQ -b $BQ $infile $ivcf > ${outpref}.counts"
	# "${MAYAROOT}/src/bin/aseQuantSample.sh -m $MQ -b $BQ -i $infile -o $outpref -v $ivcf -k $imask $CLEAN"
    else
	echo "Skipping $sample. Output already exists." 1>&2; continue;
    fi
done < "${LIST:-/proc/${$}/fd/0}"
