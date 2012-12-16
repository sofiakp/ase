#!/bin/bash

SAMPLE=SNYDER_HG19_GM12878_RNA_1_reconcile.dedup
INDIR=${MAYAROOT}/rawdata/mapped/rna/bam/personal/reconcile/dedup/q30
OUTDIR=${MAYAROOT}/rawdata/test/cuffdiff/${SAMPLE}_mat_vs_pat_v8
if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi
bsub -J ${SAMPLE}_cuffdiff -o /dev/null -eo ${OUTDIR}/cuffdiff.err -n 4 -q research-rh6 -W 96:00 -M 12288 -R "rusage[mem=12288]" "cuffdiff -o $OUTDIR -L maternal,paternal,total -p 4 --library-type fr-firststrand --no-update-check ${MAYAROOT}/rawdata/transcriptomes/gencode.v13.annotation.noM.gtf ${INDIR}/${SAMPLE}.maternal.bam ${INDIR}/${SAMPLE}.paternal.bam ${INDIR}/${SAMPLE}.bam"
#-b ${MAYAROOT}/rawdata/genomes/GM12878/GM12878.paternal.fa