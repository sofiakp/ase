#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Counts reads with the REF and ALT alleles at the positions specified.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input file.
   -o DIR [Required] Prefix of output files.
   -v VCF [Required] VCF file with positions where allele counts will be computed.
   -k FILE Mask file (R data file) with information about the variants [default: same as VCF].
   -m INT Minimum mapping quality to consider a read. [30]
   -b INT Minimum base quality to consider a position in a read. [20]
   -c     Overwrite [0]
EOF
}

CLEAN=''
INFILE=
OUTPREF=
VCF=
MASK=''
SAMPLE=''
MQ=30
BQ=20
while getopts "hi:o:v:m:b:k:s:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INFILE=$OPTARG;;
	o)
	    OUTPREF=$OPTARG;;
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

if [[ -z $OUTPREF || -z $INFILE || -z $VCF ]]; then
    usage; exit 1;
fi
if [ ! -f $INFILE ]; then
    echo 'Input file does not exist' 1>&2; exit 1;
fi
if [ ! -f $VCF ]; then
    echo 'VCF file does not exist' 1>&2; exit 1;
fi
if [[ -z $MASK ]]; then
    MASK=${VCF/.VCF/.r}
fi
if [ ! -f $MASK ]; then
    echo 'Mask file does not exist' 1>&2; exit 1;
fi

if [[ $CLEAN == "-c" || ! -f ${OUTPREF}.counts ]]; then
    ${MAYAROOT}/src/ase_cpp/bin/Ase asequantmultirg -m $MQ -b $BQ $INFILE $VCF > ${OUTPREF}.counts
fi
#if [[ $CLEAN == "-c" || ! -f ${OUTPREF}.counts.RData ]]; then
    # egrep "#CHROM" ${outpref}.counts > $tmpfile 
    # For now assume that both good and bad files are given
    # awk 'BEGIN{OFS="\t"}{if($1 != "#CHROM") print $1,$2,".",$3,$4,".",".",".",".",$5,$6,$7,$8,$9,$10,$11,$12}' ${outpref}.counts | intersectBed -a stdin -b $GOODFILE -wa | windowBed -a stdin -b $BADFILE -w 100 -v | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$5,$10,$11,$12,$13,$14,$15,$16}' >> $tmpfile
    #Rscript ${MAYAROOT}/src/rscripts/createAseData.r -i=${OUTPREF}.counts -m=$MASK -o=${OUTPREF}.counts.RData
    # rm $tmpfile
#fi