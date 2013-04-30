#!/bin/bash

usage()
{
cat <<EOF
usage: `basename $0` options
Batch signal extraction around merged peaks. Reads lines from STDIN in the form:
sname iname cell_line mark
Finds file <sname>_VS_<iname>*mat in the input directory. This is the signal track.
Then finds peak file SNYDER_HG19_<mark>_merged.encodePeak.gz in the peak directory. 
Finally, extracts the signal of the signal track in the peak regions.
OPTIONS:
   -h     Show this message and exit
   -o DIR [Required] Output directory of extracted signals.
   -i DIR [Required] Input directory of signal tracks.
   -p DIR [Required] Peak directory, or peak file. If it is a directory, then the 
          signal will be extracted around \${PDIR}/SNYDER_HG19_\${mark}_merged.bed.gz,
          otherwise, it will be extracted around the given file.
   -m STR Metafunctions. [Default: ",'mf', 'samplerate', 'mp', 10, 'ms', 150"]
   -c     Overwrite [0].
EOF
}

INDIR=
PDIR=
META=", 'mf', 'samplerate', 'mp', 10, 'ms', 150"
OUTDIR=
CLEAN=0
while getopts "hi:o:m:p:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	i)
	    INDIR=$OPTARG;;
	o)
	    OUTDIR=$OPTARG;;
	p)
	    PDIR=$OPTARG;;
	c)
	    CLEAN=1;;
	m)
	    META=$OPTARG;;
	?)
	    usage
	    exit 1;;
    esac	    
done

if [[ -z $INDIR || -z $OUTDIR || -z $PDIR ]]; then
    usage; exit 1;
fi
if [ ! -d $INDIR ]; then echo "Input directory does not exist" >&2; exit 1; fi
if [[ ! -d $PDIR && ! -f $PDIR ]]; then
    echo "PDIR does not exist" >&2; exit 1; 
fi

if [ ! -d $OUTDIR ]; then
    mkdir -p $OUTDIR
fi

JOBGRPID="/extractSignal$RANDOM"
bgadd -L 30 $JOBGRPID

while read -r sname iname cell mark; do
    # Do this if stub name is used as prefix
    sname=`echo $sname | sed -r 's/_[0-9\.]+_reconcile.dedup//'`
    #sname=`echo $sname | sed -r 's/_reconcile.dedup//'`    
    # The cell line is not used...
    # if [[ $(ls ${INDIR} | egrep ${sname}_VS_${iname}.*mat$ | wc -l) -ne 1 ]]; then
    if [[ $(ls ${INDIR} | egrep ${sname}.norm5.rawsignal.mat$ | wc -l) -ne 1 ]]; then
	echo "Input file for ${sname} not found. Skipping." >&2
	continue
    fi
    infile=${INDIR}/$(ls ${INDIR} | egrep ${sname}.norm5.rawsignal.mat$) #_VS_${iname}.*mat$)
    if [ -f $PDIR ]; then
	peakFile=$PDIR
	outpref=$(basename $PDIR)
	outpref=$(echo $outpref | sed 's/.bed[.gz]*$//')_AT_$sname
    else
	peakFile=${PDIR}/SNYDER_HG19_${mark}_merged.bed.gz
	if [[ ! -f $peakFile ]]; then 
	    echo "Peak file $peakFile not found. Skipping." >&2
	    continue
	fi
	outpref=SNYDER_HG19_${mark}_merged_AT_${sname}
    fi
    outfile=${OUTDIR}/${outpref}.mat
    if [[ $CLEAN -eq 1 || ! -f $outfile ]]; then
	script=${OUTDIR}/${outpref}.m
	echo "addpath('${SRCDIR}/signalgeneration/extractsignal/src');" > $script
	echo "addpath('${SRCDIR}/signalgeneration/extractsignal/src/fileio');" >> $script
	echo "addpath('${SRCDIR}/signalgeneration/extractsignal/src/metaFuncs');" >> $script
	echo "extractSignal('${peakFile}', '${infile}', 'if', 'bed', 'o', '${outfile}', 'ov', 'signal' $META);" >> $script
	bsub -g $JOBGRPID -J ${outpref} -o /dev/null -eo ${OUTDIR}/${outpref}.err -q research-rh6 -W 24:00 -M 12288 -R "rusage[mem=12288]" "matlab -nodisplay -nosplash < $script"
    fi
done