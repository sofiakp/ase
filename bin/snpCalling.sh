#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Call variants on a merged set of BAM files.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input list of fles.
   -o DIR [Required] Output directory for merged BAMs.
   -v DIR [Required] Output directory for variant calls. 
   -s STR [Required] Sample.
   -r STR [Rqquired] Chromosome (or region).
   -f FILE Fasta file.
   -c     Overwrite [0]
EOF
}

CLEAN=0
INDIR=
OUTDIR=
VARDIR=
SAMPLE=
CHROM=
FASTA=${MAYAROOT}/rawdata/genomes/encodeHg19Male/encodeHg19Male.fa
while getopts "hi:o:v:s:r:f:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	v) 
	    VARDIR=$OPTARG;;
	s)
	    SAMPLE=$OPTARG;;
	r)
	    CHROM=$OPTARG;;
	f)
	    FASTA=$OPTARG;;
	c) 
	    CLEAN=1;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || -z $VARDIR || -z $OUTDIR || -z $SAMPLE || -z $CHROM ]]; then
    usage; exit 1;
fi
if [ ! -f $FASTA ]; then
    echo "Input sequence file does not exist" 1>&2; exit 1;
fi
if [[ ! -d $OUTDIR ]]; then
    mkdir -p $OUTDIR
fi
if [[ ! -d $VARDIR ]]; then
    mkdir -p $VARDIR
fi

files=$INDIR

outpref=${OUTDIR}/${SAMPLE}_${CHROM}
varpref=${VARDIR}/${SAMPLE}_${CHROM}

if [[ $CLEAN -eq 1 || ! -f ${outpref}.bam.bai ]]; then
    echo $files > ${outpref}_files.txt
    $SAMTOOLS18 merge -R $CHROM ${outpref}_rg.bam $files
    java -Xmx8g -jar ${PICARD}/AddOrReplaceReadGroups.jar I=${outpref}_rg.bam O=${outpref}.bam RGPL=Illumina RGSM=$SAMPLE RGLB=$SAMPLE RGPU=PU VALIDATION_STRINGENCY=LENIENT
    rm ${outpref}_rg.bam
    $SAMTOOLS18 index ${outpref}.bam
fi

if [[ $CLEAN -eq 1 || ! -f ${varpref}_raw.vcf ]]; then
    $SAMTOOLS18 mpileup -f $FASTA -q 20 -Q 15 -u ${outpref}.bam | $BCFTOOLS18 view -vcg - > ${varpref}_raw.vcf
fi
