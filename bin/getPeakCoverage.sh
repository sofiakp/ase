#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Computes coverage (i.e. # peaks) of ChIP datasets (all replicates + merged) on the corresponing peaks.
OPTIONS:
   -h     Show this message and exit
   -i DIR [Required] Directory with BED files.
   -p DIR [Required] Directory with peak files OR bed file. In the second case, coverage will be computed on the provided regions.
   -o DIR [Required] Output directory.
   -c     Overwrite [0].
   -m     Compute coverage on the merged peaks. By default coverage is computed on individual specific 
          peaks and then the counts for all individul-specific peak overlapping the same merged peak are summed up.
EOF
}

INDIR=
PDIR=
OUTDIR=
CLEAN=0
MERGED=0
while getopts "hi:o:p:cm" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	p)
	    PDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	c)
	    CLEAN=1;;
	m)
	    MERGED=1;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $PDIR || -z $INDIR || -z $OUTDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $INDIR ]; then echo "Input directory does not exist" >&2; exit 1; fi
if [[ ! -d $OUTDIR ]]; then mkdir -p $OUTDIR; fi

JOBGRPID="/getPeakCoverage$RANDOM"
bgadd -L 30 $JOBGRPID 

while read inputline; do
    chip=$(echo ${inputline} | awk '{print $1}' | sed -r 's/\.bam/\.bedpe.gz/g') # get names of BEDPE files
    shortstub=$(echo ${inputline} | awk '{print $3}' | awk 'BEGIN{FS="_"}{print $1"_"$2"_"$4}' ) # Prefix of merged peaks, eg. SNYDER_HG19_H3K4ME3

    # Get the peak file corresponding to that input line
    if echo ${chip} | grep -q ';'; then
	chipstub=$(echo ${inputline} | awk '{print $3}') # If multiple chip files then use column 3 in file as ChIPstub
	#mergechip=
    else
	chipstub=$(echo ${chip} | sed -r 's/\.bedpe.gz//g') # If one chip file then use remove .bed.gz from the file name and use that as ChIPstub
	#mergechip=${INDIR}/${chipstub}_0_reconcile.dedup.bedpe.gz
    fi
    controlstub=$(echo ${inputline} | awk '{print $2}' | sed -r -e 's/;.*$//g' -e 's/\.bam//g') # If multiple inputs, first one is used as stub

    # For each replicate in the first column of the input line, i.e. for each replicate that was used to
    # get the peakfile, compute the coverage of the replicate in the peaks.
    # There are three alternatives here:
    # A) If PDIR is not a file and MERGED is 1: The coverage is computed on the peaks of that input line (i.e. individual-specific peaks) a
    # nd then counts for individual-specific peaks are aggregated over all such peaks overlapping the same "merged" peak.
    # B) If PDIR is not a file a MERGED is 0: The coverage is compute on the merged peaks.
    # C) If PDIR is a file: The coverage is computed on the regions of this file.
    for c in $(echo $chip | sed 's/;/\n/g'); do
	bedfile=${INDIR}/${c}
	if [ ! -f $bedfile ]; then echo "Missing file $bedfile" >&2; continue; fi
	if [ -f $PDIR ]; then
	    mergepeaks=$PDIR
	    outpref=$(basename $PDIR)
	    outpref=${outpref/.bed}_AT_${c/\.bedpe\.gz}.bed
	else
	    mergepeaks=${PDIR}/${shortstub}_merged.bed.gz
	    if [ ! -f $mergepeaks ]; then echo "Missing file $mergepeaks" >&2; continue; fi
	    
	    outpref=${shortstub}_merged_AT_${c/\.bedpe\.gz/}.bed
	fi
	outfile=${OUTDIR}/$outpref
	#script=${OUTDIR}/${outpref}.sh
	#echo "#!/bin/bash" > $script
	#todo=0

	if [[ $CLEAN -eq 1 || ! -f $outfile ]]; then
	    if [[ -f $PDIR || $MERGED -eq 1 ]]; then
		bsub -q research-rh6 -g "${JOBGRPID}" -J "$outpref" -W 24:00 -o /dev/null -e /dev/null "coverageBed -a $bedfile -b $mergepeaks -counts | sort -V > $outfile" 
	    else
		peakFile=${PDIR}/${chipstub}_VS_${controlstub}_peaks.encodePeak.gz
		if [ ! -f $peakFile ]; then echo "Missing file $peakFile" >&2; continue; fi
		# Get the coverage of reads on the peaks. Then overlap these peaks with the merged set of peaks (over all individuals)
	        # and aggregate the counts for all individual-specific peaks overlapping each merged peak.
		# bsub -q research-rh6 -g "${JOBGRPID}" -J "$outpref" -W 24:00 -o /dev/null -e /dev/null "coverageBed -a $bedfile -b $peakFile -counts | intersectBed -a $mergepeaks -b stdin -wa -wb -loj | python $MAYAROOT/src/python/sumUpCounts.py > $outfile" 
	    fi
#	else
#	    echo "$outpref exists. Skipping." >&2
	fi
    done
    # TODO: check if merged BED file exists and compute merged coverage
done
