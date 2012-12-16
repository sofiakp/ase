#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Concatenate stat files from MarkDuplicates.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -r     Does not have read groups (maternal, paternal, ambiguous).
   -s     Single end. [0]
EOF
}
INDIR=
RG=1
SINGLE=0
while getopts "hi:sr" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	r)
	    RG=0;;
	s)
	    SINGLE=1;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $INDIR ]; then
    echo 'Indir does not exist' 1>&2; exit 1;
fi

if [[ $SINGLE -eq 1 ]]; then
    col=2
else
    col=3
fi
if [[ $RG -eq 1 ]]; then
    echo -e "File\tPat_pairs\tPat_dup_rate\tMat_pairs\tMat_dup_rate\tAmb_pairs\tAmb_dup_rate"
    for i in `ls $INDIR | egrep dedup.stats`; do
	base=`basename $i`
	base=${base/.stats/}
	p=`egrep "^paternal" ${INDIR}/$i | cut -f2,3,4,8 | awk 'BEGIN{OFS="\t"}{print $1 + 2*$2, $4}'`
	m=`egrep "^maternal" ${INDIR}/$i | cut -f2,3,4,8 | awk 'BEGIN{OFS="\t"}{print $1 + 2*$2, $4}'`
	a=`egrep "^ambiguous" ${INDIR}/$i | cut -f2,3,4,8 | awk 'BEGIN{OFS="\t"}{print $1 + 2*$2, $4}'`
	echo -e "${base}\t${p}\t${m}\t${a}"
    done
else
    echo -e "File\tPairs\tDup_rate"
    for i in `ls $INDIR | egrep dedup.stats`; do
	base=`basename $i`
	base=${base/.stats/}
	p=`head -8 ${INDIR}/$i | tail -1 | cut -f2,3,4,8 | awk 'BEGIN{OFS="\t"}{print $1 + 2*$2, $4}'`
	echo -e "${base}\t${p}"
    done
fi