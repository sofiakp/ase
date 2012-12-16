#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Call variants on a set of individuals. All bam files for each individual are merged and variants are called on all individuals at the same time (multisample calls).
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Input directory with bam files.
   -o DIR [Required] Output directory for merged BAMs. All BAMs for each individual will be merged and then split by chromosome.
   -v STR [Required] Prefix for output files. Output files will be <prefix>_chr*_raw.vcf.
   -s STR [Required] Prefixes of BAM files for each individual, separated by ;. For eg. "SNYDER_HG19_foo;SNYDER_HG19_bar": in this case all BAM files in the input directory starting with SNYDER_HG19_foo will be merged together and then split by chromosome. The same will be done for the second prefix provided. Then variants will be called on all merged files for the same chromosome.
   -f FILE Fasta file [Default \$MAYAROOT/rawdata/genomes/encodeHg19Male/encodeHg19Male.fa].
   -c     Overwrite [0]
EOF
}

CLEAN=''
INDIR=
OUTDIR=
VARPREF=
INDIV=
FASTA=${MAYAROOT}/rawdata/genomes/encodeHg19Male/encodeHg19Male.fa
while getopts "hi:o:v:s:f:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	v) 
	    VARPREF=$OPTARG;;
	s)
	    INDIV=$OPTARG;;
	f)
	    FASTA=$OPTARG;;
	c) 
	    CLEAN='-c';;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || -z $VARPREF || -z $OUTDIR || -z $INDIV ]]; then
    usage; exit 1;
fi
if [ ! -d $INDIR ]; then
    echo "Input directory does not exist" 1>&2; exit 1;
fi
if [ ! -f $FASTA ]; then
    echo "Input sequence file does not exist" 1>&2; exit 1;
fi
if [[ ! -d $OUTDIR ]]; then
    mkdir -p $OUTDIR
fi

for c in {1..22} X Y; do 
    script=${VARPREF}_chr${c}.sh
    echo "#!/bin/bash" > $script
    
    # Merge all BAM files in INDIR with the same prefix (eg. BAM files for different marks).
    # Then split the merged files by chromosome for parallel variant calling.
    for s in $(echo $INDIV | sed -r 's/;/\n/g'); do
	files=$(find $INDIR -maxdepth 1 -name ${s}*bam -printf "%p ")
	merged=${OUTDIR}/${s}_chr${c}
	if [[ $CLEAN -eq 1 || ! -f ${merged}.bam ]]; then
	    echo "$SAMTOOLS18 merge -R chr$c ${merged}_rg.bam $files" >> $script
	    echo "java -Xmx8g -jar ${PICARD}/AddOrReplaceReadGroups.jar I=${merged}_rg.bam O=${merged}.bam RGPL=Illumina RGSM=$s RGLB=$s RGPU=PU VALIDATION_STRINGENCY=LENIENT" >> $script
	    echo "rm ${merged}_rg.bam" >> $script
	    echo "$SAMTOOLS18 index ${merged}.bam" >> $script
	fi
    done

    # Call variants on all individuals.
    outfile=${VARPREF}_chr${c}_raw.vcf
    if [[ $CLEAN -eq 1 || ! -f $outfile ]]; then
	infiles=$(echo $INDIV | sed -r 's/;/\n/g' | xargs -I fname echo ${OUTDIR}/fname_chr${c}.bam | sed -e :a -e '/$/N;s/\n/ /;ta')
	echo "$SAMTOOLS18 mpileup -f $FASTA -q 20 -Q 15 -u $infiles | $BCFTOOLS18 view -vcg - > ${outfile}" >> $script
	bsub -q research-rh6 -J chr${c}_snps -W 96:00 -M 12288 -R "rusage[mem=12288]" -o /dev/null -e /dev/null < ${script}
	#rm $script
    fi
done
