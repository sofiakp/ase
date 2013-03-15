#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Aligns an RNA sample with TopHat
OPTIONS:
   -h            Show this message and exit
   --fq1 FILE    [Required]
   --fq2 FILE    [Required] Fastq files
   --bamdir DIR  [Required] Dir where bam files will be written
   --mpref STR   [Required] Prefix of name sorted file.
   --seqpref STR [Required] Prefix for genome index files
   --trpref STR  [Required] Prefix for transcriptome index files (will be created if they don't exist)
   --gtf FILE    [Required] GTF file with known transcriptome.
   -c            Overwrite output files [0]
EOF
}

ARGS=`getopt -o "hc" -l "fq1:,fq2:,bamdir:,seqpref:,trpref:,mpref:,gtf:" -- "$@"`
eval set -- "$ARGS"

CLEAN=0
BAMDIR=
FQ1=
FQ2=
SEQPREF=
TRPREF=
MPREF=
GTF=
while [ $# -gt 0 ]; do
    case $1 in
	-h) usage; exit;;
	--fq1) FQ1=$2; shift 2;;
	--fq2) FQ2=$2; shift 2;;
	--bamdir) BAMDIR=$2; shift 2;;
	--seqpref) SEQPREF=$2; shift 2;;
	--mpref) MPREF=$2; shift 2;;
	--trpref) TRPREF=$2; shift 2;;
	--gtf) GTF=$2; shift 2;;
	-c) CLEAN=1; shift;;
	--) shift; break;;
    esac	    
done

if [ $# -ne 0 ]; then
    usage; exit 1;
fi

if [[ -z $BAMDIR || -z $FQ1 || -z $FQ2 || -z $SEQPREF || -z $TRPREF || -z $MPREF || -z $GTF ]]; then
    usage; exit 1;
fi

if [ ! -d $BAMDIR ]; then
    mkdir -p $BAMDIR
fi

# TopHat options
ALN="-p 4 --no-discordant --library-type fr-firststrand --b2-sensitive" #--no-novel-juncs"

if [[ -s ${FQ1} && -s ${FQ2} ]]; then
    # Check Illumina version
    format=`zcat $FQ1 | python ${MAYAROOT}/src/python/checkFastqFormat.py`
    if [[ $format == '1' ]]; then
	ALN="$ALN --solexa1.3-quals"
    fi
    if [[ -f ${TRPREF}.gff ]]; then
	if [[ ! -f ${TRPREF}.1.bt2 ]]; then
	    # TopHat will break if the gff is there but the index is not
	    rm ${TRPREF}.gff
	    ALN="$ALN --GTF $GTF"
	fi
    else
	ALN="$ALN --GTF $GTF"
    fi
    zcat $FQ1 | sed -r 's/_[12]:N:0:[ACGT]+//' | gzip -c > ${FQ1}_TMP.gz
    zcat $FQ2 | sed -r 's/_[12]:N:0:[ACGT]+//' | gzip -c > ${FQ2}_TMP.gz
    if [[ $CLEAN -eq 1 || ! -f ${BAMDIR}/accepted_hits.bam ]]; then
	tophat2 $ALN -o $BAMDIR -r 200 --mate-std-dev 20 --transcriptome-index $TRPREF $SEQPREF ${FQ1}_TMP.gz ${FQ2}_TMP.gz
    fi

    if [[ $CLEAN -eq 1 || ! -f ${MPREF}.bam ]]; then
	$SAMTOOLS18 view -H ${BAMDIR}/accepted_hits.bam | sed 's/SO:coordinate/SO:unsorted/' > ${BAMDIR}/header.sam
	num=`$SAMTOOLS18 view ${BAMDIR}/accepted_hits.bam | head -1 | egrep "_[12]:N:0:[ACGT]+" | wc -l`
	if [[ $num -gt 0 ]]; then
	    $SAMTOOLS18 cat -h ${BAMDIR}/header.sam ${BAMDIR}/accepted_hits.bam ${BAMDIR}/unmapped.bam | $SAMTOOLS18 view -h -F0x100 - | sed -r 's/_[12]:N:0:[ACGT]+//' | $SAMTOOLS18 view -Sb - | $SAMTOOLS18 sort -n -m 2000000000 - $MPREF
	else
	    $SAMTOOLS18 cat -h ${BAMDIR}/header.sam ${BAMDIR}/accepted_hits.bam ${BAMDIR}/unmapped.bam | $SAMTOOLS18 view -b -F0x100 - | $SAMTOOLS18 sort -n -m 2000000000 - $MPREF
	fi
    fi
fi
