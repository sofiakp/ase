#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Counts reads on exons.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input file.
   -o DIR [Required] Prefix of output files.
   -g FILE [Required] GTF file with gene annotations.
   -e DIR Prefix of exon count files. If missing, exon counting won't be run.
   -m INT Minimum mapping quality to consider a read. [30]
   -c     Overwrite [0]
EOF
}

CLEAN=''
INFILE=
OUTPREF=
EXPREF=
GTF=
DIRECTION="no"
SINGLE="-p yes"
MQ=30
while getopts "hi:o:e:g:m:d:s:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INFILE=$OPTARG;;
	o)
	    OUTPREF=$OPTARG;;
	e)
	    EXPREF=$OPTARG;;
	g) 
	    GTF=$OPTARG;;
	m)
	    MQ=$OPTARG;;
	s)
	    SINGLE="-p no";;
	d) 
	    DIRECTION=$OPTARG;;
	c) 
	    CLEAN='-c';;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $OUTPREF || -z $INFILE || -z $GTF ]]; then
    usage; exit 1;
fi
if [ ! -f $INFILE ]; then
    echo 'Input file does not exist' 1>&2; exit 1;
fi
if [ ! -f $GTF ]; then
    echo 'GTF file does not exist' 1>&2; exit 1;
fi

GFF=${GTF/.gtf/.flat.gff}
if [[ ! -z $EXPREF && ! -f $GFF ]]; then
    echo 'GFF file does not exist' 1>&2; exit 1;
fi  

if [[ $CLEAN == "-c" || ! -f ${OUTPREF}_ambiguous.genecounts ]]; then
    for g in ambiguous paternal maternal; do
	# Input file MUST be sorted by name
	# Select reads that are mapped and have a mapped mate.
	# Tophat doesn't set the "properly paired" bit correctly, so selecting for that would remove too many pairs.
	samtools view -F0xC -r $g $INFILE | htseq-count -m intersection-nonempty -s $DIRECTION -a $MQ -t exon -i gene_id - $GTF > ${OUTPREF}_${g}.genecounts
    done
fi

#if [[ ! -z $EXPREF && ($CLEAN == "-c" || ! -f ${EXPREF}_ambiguous.exoncounts) ]]; then
#    for g in ambiguous paternal maternal; do
#	samtools view -F0xC -r $g $INFILE | python $MAYAROOT/tools/DEXSeq/inst/python_scripts/dexseq_count.py -s $DIRECTION -a $MQ $SINGLE $GFF - ${EXPREF}_${g}.exoncounts
#    done
#fi
