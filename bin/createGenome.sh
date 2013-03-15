#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Creates personal genome.
OPTIONS:
   -h     Show this message and exit
   -v VCF [Required] VCF file
   -i STR [Required] Individual name (as in VCF)
   -o DIR Output directory. Default is \$MAYAROOT/rawdata/genomes/<indiv>
   -d STR Genome dictionary. Default is \$MAYAROOT/rawdata/genomes/encodeHg19Male/hg19noM.fa.fai.
   -s DIR Dir with fa files, one per chromosome. Default is \$SEQDIR/encodeHg19Male
   -f     Set this option if the individual is a female.
   -r     Set this option if the VCF does NOT have chr in the chromosome names.
EOF
}

VCF=
INDIV=
OUTDIR=
DICT=${MAYAROOT}/rawdata/genomes/encodeHg19Male/hg19noM.fa.fai
SDIR=${SEQDIR}/encodeHg19Male
FEMALE=''
CHROM=''
while getopts "hi:o:d:v:fr" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIV=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	d)
	    DICT=$OPTARG;;
	v)
	    VCF=$OPTARG;;
	f) 
	    FEMALE='-f';;
	r)
	    CHROM='-c';;
	s) 
	    SDIR=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $VCF || -z $INDIV ]]; then
    usage; exit 1;
fi
if [[ -z $OUTDIR ]]; then
    OUTDIR=${MAYAROOT}/rawdata/genomes/$INDIV
fi
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi
if [ ! -f $DICT ]; then
    echo 'Genome dictionary does not exist' 1>&2; exit 1;
fi
if [ ! -d $SDIR ]; then
    echo 'Genome dictionary does not exist' 1>&2; exit 1;
fi

if [[ $INDIV =~ NA[0-9]+ ]]; then
    outindiv=${INDIV/NA/GM}
else
    outindiv=$INDIV
fi

if [[ $INDIV == "SNYDER" ]]; then
    #bsub -J $outindiv -e /dev/null -o /dev/null -n 1 -q research-rh6 -W 2:00 -M 8192 -R "rusage[mem=8192]" 
    cat $VCF | python ${MAYAROOT}/src/python/addSnyderSnps.py $SDIR $DICT ${OUTDIR}/${outindiv}
else
    bsub -J $outindiv -eo ${OUTDIR}/${outindiv}.err -o /dev/null -n 1 -q research-rh6 -W 2:00 -M 8192 -R "rusage[mem=8192]" "python ${MAYAROOT}/src/python/addSeqphaseSnpsToFa.py $SDIR $DICT ${MAYAROOT}/rawdata/phasing/San/merged/CG_AFR_ ${OUTDIR}/${outindiv} $INDIV --seqphase $VCF $FEMALE"
    #bsub -J $outindiv -e /dev/null -o /dev/null -n 1 -q research-rh6 -W 2:00 -M 8192 -R "rusage[mem=8192]" "cat $VCF | python ${MAYAROOT}/src/python/addSnpsToFa.py $SDIR $DICT ${OUTDIR}/${outindiv} $INDIV --unphased ${OUTDIR}/${outindiv}.unphased.txt $FEMALE $CHROM"
fi