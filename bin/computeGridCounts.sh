#!/bin/bash

INDIR=$1
OUTDIR=$2

while read -r sample; do
    sample=${sample^^} # This will convert to uppercase
    if [[ "$sample" =~ ^[0-9]+ ]]; then # Correct HapMap names
	sample="GM"$sample
    fi
    sample="SNYDER_HG19_${sample}"
    infile=${INDIR}/${sample}_reconcile.dedup.bam
    outfile=${OUTDIR}/${sample}_grid_chr15.bed
    bsub -J "${sample}_grid" -e /dev/null -o /dev/null "samtools view -b $infile chr15 | coverageBed -abam stdin -b $MAYAROOT/rawdata/genomes/chr15_grid.bed -counts | sort -V > $outfile"
done
