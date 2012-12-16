#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Run cuffdiff on paternal vs maternal samples.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory.
   -o DIR [Required] Output directory.
   -s STR [Required] Individual.
   -m STR [Required] Mark. 
   -g FILE GTF file with gene annotations [Default: \$MAYAROOT/rawdata/transcriptomes/gencode.v13.annotation.gtf]
   -f FILE Fasta sequence file [Default: \${MAYAROOT}/rawdata/genomes/hg19/encodeHg19Male.fa]
   -c     Overwrite output directory.
EOF
}

INDIR=
OUTDIR=
INDIV=
MARK=
CLEAN=0
GTF=${MAYAROOT}/rawdata/transcriptomes/gencode.v13.annotation.gtf
FASTA=${MAYAROOT}/rawdata/genomes/hg19/encodeHg19Male.fa
while getopts "hi:o:s:m:g:f:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	s)
	    INDIV=$OPTARG;;
	m) 
	    MARK=$OPTARG;;
	g)
	    GTF=$OPTARG;;
	f)
	    FASTA=$OPTARG;;
	c) 
	    CLEAN=1;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || -z $OUTDIR || -z $INDIV || -z $MARK ]]; then
    usage; exit 1;
fi

if [ ! -f $GTF ]; then
    echo "GTF file does not exist." 1>&2; exit 1;
fi
if [ ! -f $FASTA ]; then
    echo "FASTA file does not exist." 1>&2; exit 1;
fi
if [ ! -d $INDIR ]; then
    echo "Input directory does not exist." 1>&2; exit 1;
fi

if [[ "$INDIV" =~ ^[0-9]+$ ]]; then # Correct HapMap names
    INDIV="GM"$INDIV
fi

dirpref=SNYDER_HG19_${INDIV}_${MARK}
cuffdir=${OUTDIR}/${dirpref}_paternal_vs_maternal

pats=
for d in `ls $INDIR | egrep ${dirpref}.*_paternal`; do
    bamfile=${INDIR}/${d}/accepted_hits.bam
    if [ -f $bamfile ]; then
	echo $d
	if [ -z $pats ]; then
	    pats=$bamfile
	else
	    pats="${pats},$bamfile"
	fi
    fi
done

mats=
for d in `ls $INDIR | egrep ${dirpref}.*_maternal`; do
    bamfile=${INDIR}/${d}/accepted_hits.bam
    if [ -f $bamfile ]; then
	echo $d
	if [ -z $mats ]; then
	    mats=$bamfile
	else
	    mats="${mats},$bamfile"
	fi
    fi
done

if [[ ! -z $pats && ! -z $mats ]]; then
    if [ ! -d $cuffdir ]; then 
	mkdir -p $cuffdir
    elif [[ $CLEAN -eq 0 && -f ${cuffdir}/run.info ]]; then
	echo "Output directory exists. Skipping." 1&>2; exit
    fi
    labels="paternal,maternal"
    FASTA=${MAYAROOT}/rawdata/genomes/${INDIV}/${INDIV}.paternal.fa
    bsub -J ${dirpref}_cuffdiff -o /dev/null -eo ${cuffdir}/cuffdiff.err -n 6 -q research-rh6 -W 96:00 -M 12288 -R "rusage[mem=12288]" "cuffdiff -o $cuffdir -L $labels -p 6 -N -u --library-type fr-firststrand --no-update-check $GTF $pats $mats"
    # -b $FASTA
else 
    echo "Paternal or maternal files (or both) are missing." 1>&2; exit 1;
fi