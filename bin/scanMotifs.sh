#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Runs BWA on a single sample
OPTIONS:
   -h            Show this message and exit
   -o Output directory
   -m FILE motif file
   -i STR [Optional] motif id 
   -c Overwrite
EOF
}

CLEAN=0
ODIR=
MOTIFS=
MID=
while getopts "ho:m:i:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	o)
	    ODIR=$OPTARG;;
	c) 
	    CLEAN=1;;
	m)
	    MOTIFS=$OPTARG;;
	i)
	    MID=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $ODIR || -z $MOTIFS ]]; then
    usage; exit 1;
fi

if [ ! -d $ODIR ]; then
    mkdir -p $ODIR
fi

if [[ -z $MID ]]; then
    SELMOT=""
else
    SELMOT="--motif $MID"
fi

while read indiv; do
    indiv=${indiv^^}
    if [[ "$indiv" =~ ^[0-9]+$ ]]; then # Correct HapMap names
	indiv="GM"$indiv
    fi
    seqpref
	if [[ ! -f ${seqpref}.maternal.fa || ! -f ${seqpref}.paternal.fa ]]; then
	    seqpref=${seqpref}/${indiv}
	    echo $seqpref
	    if [[ ! -f ${seqpref}.maternal.fa || ! -f ${seqpref}.paternal.fa ]]; then
		echo "Could not file genome files for $sample. Skipping..." 1>&2 ; continue;
	    fi
	fi
    for p in maternal paternal; do
    bsub -J fimo_test -eo $outdir/job.err -o /dev/null -W 96:00 -M 12288 -R "rusage[mem=12288]" "fimo --bgfile motif-file --oc $outdir --max-stored-scores 1000000 --output-pthresh 0.0001 CTCF.meme $MAYAROOT/rawdata/genomes/GM12878/GM12878.maternal.fa"
    done
done < "${LIST:-/proc/${$}/fd/0}"
