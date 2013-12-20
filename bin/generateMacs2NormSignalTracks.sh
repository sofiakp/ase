#!/bin/bash

usage()
{
cat <<EOF
Usage: $(basename $0) [options] <iDir> <oDir>
<iDir>: input directory whose subdirectories contain mapped ChIP and control data.
<oDir>: output directory. Log files will go to <oDir>/logFiles and bedgraph files will go to <oDir>/temp. 

OPTIONS:
   -h      Show this message and exit
   -p FILE Input file containing pairs of ChIP (column 1), control dataset names (column 2), and a composite experiment name (column3). 
           Multiple ChIP and control files can be present seperated by ; in Col1 and Col2. They will be merged before passing them to MACS.
   -g NUM  Genome size. Default: hs.
   -s FILE Tab delimited chrName\tsize: Default : \$GENOMESIZEDIR/hg19.genome
   -m NUM  Memory limit, if non-zero then it represents X GB. Set to 0 to use the default.
   -f FILE File containing estimated fragment length for each dataset. Col1: dataset name [tab] Col2:frag length
   -c      Overwrite.
EOF
}

IFILE=
GENOMESIZE='hs'
MEM=0
CHRSIZE="${GENOMESIZEDIR}/hg19.genome"
FRAGLENFILE=
CLEAN=0
while getopts "hp:g:s:m:f:c" opt
do
    case $opt in
	h)
	    usage; exit;;
	p)
	    IFILE=$OPTARG;;
	g)
	    GENOMESIZE=$OPTARG;;
	s)
	    CHRSIZE=$OPTARG;;
	m)
	    MEM=$OPTARG;;
	f)
	    FRAGLENFILE=$OPTARG;;
	c) 
	    CLEAN=1;;
	?)
	    usage
	    exit 1;;
    esac	    
done
shift $(( OPTIND-1 ))

if [[ "$#" -ne 2 ]]; then usage; exit 1; fi

IDIR=$(readlink -f $1)
if [[ ! -d "${IDIR}" ]]; then echo "ERROR: Argument <iDir> ${IDIR} does not exist" >&2 ; exit 1; fi
ODIR=$2
if [[ ! -d "${ODIR}/logFiles" ]]; then mkdir -p ${ODIR}/logFiles; fi
if [[ ! -d "${ODIR}/temp" ]]; then mkdir -p ${ODIR}/temp; fi

if [[ ! -z $IFILE && ! -f "${IFILE}" ]]; then echo "ERROR: Pair file ${IFILE} does not exist" >&2 ; exit 1; fi

if [[ ! -f "${CHRSIZE}" ]]; then echo "ERROR: File with chromosome sizes ${CHRSIZE} does not exist" >&2 ; exit 1; fi

memLimit=$(( MEM * 1024 ))

if [[ ! -z ${FRAGLENFILE} && ! -f "${FRAGLENFILE}" ]]
then
    echo "ERROR: Fragment sizes file ${FRAGLENFILE} does not exist" >&2
    exit 1
fi

# Check if you can call Rscript
if [[ -z $(which macs2) ]]; then echo 'ERROR: MACS executable not in $PATH' >&2; exit 1; fi

# Maximum number of jobs to run at a time
JOBGRPID="/histVarSignal${RANDOM}"
bgadd -L 30 $JOBGRPID

# ========================================
# Read pairFile line by line
# Create shell script
# Submit to cluster
# ========================================

while read inputline ; do
    #echo $inputline
    # extract first column as ChIP file name, replace .bam with .bed.gz and sort filenames
    chip=$(echo ${inputline} | awk '{print $1}' | sed -r 's/\.bam/.bed.gz/g' | sed 's/;/\n/g' | sort | sed -e :a -e '/$/N;s/\n/;/;ta')
    #echo $chip
    if echo ${chip} | grep -q ';'; then
	chipstub=$(echo ${inputline} | awk '{print $3}') # If multiple chip files then use column 3 in file as ChIPstub
    else
	chipstub=$(echo ${chip} | sed -r 's/\.bed\.gz//g') # If one chip file then use remove .bed.gz from the file name and use that as ChIPstub
    fi
 
    control=$(echo ${inputline} | awk '{print $2}' | sed -r 's/\.bam/.bed.gz/g') # extract second column as control file name, replace .bam with .bed.gz
    if echo ${control} | grep -q ';'; then
	controlstub='SNYDER_HG19_all_INPUT'
    else
	controlstub=$(echo ${control} | sed -r -e 's/;.*$//g' -e 's/\.bed\.gz//g') # Use first file name if multiple separated by ; and remove the .bed.gz extension to generate controlstub
    fi

    outFile="${ODIR}/logFiles/${chipstub}_VS_${controlstub}.out"
    errFile="${ODIR}/logFiles/${chipstub}_VS_${controlstub}.err"
    peakFile="${ODIR}/temp/${chipstub}_VS_${controlstub}"
    fcFile="${peakFile}.fc.signal.bedgraph.gz"
    llrFile="${peakFile}.llr.signal.bedgraph.gz"

    if [[ -e $outFile && -e $errFile && $CLEAN -eq 0 ]]; then
	if grep -q "Finished!.*llr.bedgraph" $errFile; then
	    oldChip=$(egrep "^Combining ChIP replicates" $outFile | sed 's/Combining ChIP replicates: //' | sed 's/ /\n/g' | sort | xargs -I fname basename fname | sed -e :a -e '/$/N;s/\n/;/;ta')
	    # If the basenames of the replicate files haven't changed, skip.
	    if [[ "$chip" == "$oldChip" ]]; then
		echo "Skipping ${peakFile}." >&2
		continue
	    else
		echo "Replicates for ${peakFile} have changed. Rerunning." >&2
	    fi
	else
	    echo "Log files for ${peakFile} exist, but it looks like the jobs didn't complete successfully. Rerunning." >&2
	fi
    else
	echo "Starting ${peakFile}." >&2
    fi
    #continue # For debugging

    nchip=$(echo ${chip} | sed -r 's/;/\n/g' | wc -l) # number of chip files
    chip=$(echo ${chip} | sed -r 's/;/\n/g' | sort | xargs -I fname find "${IDIR}" -name fname -printf "%p ") # separate file names with space

    # Search for full file path in IDIR
    nchipfound=$(echo ${chip} | sed -r 's/ /\n/g' | wc -l) # number of chip files that were found in IDIR
    if [[ ${nchipfound} -ne ${nchip} ]]; then echo "ERROR: Some of the ChIP files ${chip} were not found in $IDIR" >&2 ; continue ; fi
    
    ncontrol=$(echo ${control} | sed -r 's/;/\n/g' | wc -l) # number of control files
    control=$(echo "${control}" | sed -r 's/;/\n/g' | xargs -I fname find "${IDIR}" -name fname -printf "%p ") # separate file names with space
    ncontrolfound=$(echo ${control} | sed -r 's/ /\n/g' | wc -l) # number of control files that were found in IDIR
    if [[ ${ncontrolfound} -ne ${ncontrol} ]]; then echo "ERROR: Some of the control files ${control} were not found in $IDIR" >&2 ; continue ; fi
 
    # Get fragLens corresponding to ChIP files
    fraglen=0
    fcount=0
    if [[ ! -z ${FRAGLENFILE} ]]; then
      for currFile in $(echo ${chip}); do
	  currBase=$(echo $(basename ${currFile} | sed -r 's/\.bed\.gz/\\./g'))
	  if grep -q ${currBase} ${FRAGLENFILE}; then
	      currFragLen=$(grep ${currBase} ${FRAGLENFILE} | awk '{print $2}')
	      fraglen=$((fraglen + currFragLen))
	      fcount=$((fcount + 1))
	  else
	      echo "WARNING: No fragment length found corresponding to file ${currFile}" >&2
	  fi
      done
      [[ ${fcount} -eq 0 ]] && fcount=1
      fraglen=$((fraglen / (2 * fcount))) # Divide by 2 since MACS expects shift-size which is 1/2 fragment length
      if [[ ${fraglen} -eq 0 ]]; then
	  echo "ERROR: Fragment length is 0 due to missing file names in ${FRAGLENFILE}" >&2
	  continue
      fi
    fi

    scriptName="temp${RANDOM}${RANDOM}.sh" # script name
    echo '#!/bin/bash' > ${scriptName}
    echo 'tmpdir="${TMP}/tmp${RANDOM}_${RANDOM}"' >> ${scriptName}
    echo 'mkdir ${tmpdir}' >> ${scriptName}
  
    # -------------------------
    # Create temp copies of ChIP and control
    # gunzip and concatenate multiple files if necessary
    # -------------------------	
    echo 'combchip="${tmpdir}/'"${chipstub}"'"' >> ${scriptName}
    echo 'if [[ -f "${combchip}" ]]; then rm -rf "${combchip}"; fi' >> ${scriptName}
    echo "echo Combining ChIP replicates: ${chip}" >> ${scriptName}
    echo "zcat ${chip} | awk 'BEGIN{OFS="'"\t"}{$4="N";$5="1000";print $0}'"'"' >> "${combchip}"' >> ${scriptName}
    
    echo 'combcontrol="${tmpdir}/'"${controlstub}"'"' >> ${scriptName}
    echo 'if [[ -f "${combcontrol}" ]]; then rm -rf "${combcontrol}"; fi' >> ${scriptName}
    echo "echo Combining Control replicates: ${control}" >> ${scriptName}
    echo "zcat ${control} | awk 'BEGIN{OFS="'"\t"}{$4="N";$5="1000";print $0}'"'"' >> "${combcontrol}"' >> ${scriptName}
    
    # -------------------------
    # Complete script
    # -------------------------
    if [[ -z ${FRAGLENFILE} || ${fraglen} -eq 0 ]]; then
        #echo 'macs2 callpeak -t "${combchip}" -c "${combcontrol}" -f BED'" -n ${peakFile} -g ${GENOMESIZE} -p 1e-2 -m 5,30 -B --SPMR" >> ${scriptName}
	echo 'macs2 callpeak -t "${combchip}" -c "${combcontrol}" -f BED'" -n ${peakFile} -g ${GENOMESIZE} -p 1e-2 --nomodel --shiftsize 73 -B --SPMR" >> ${scriptName}
    else
	echo 'macs2 callpeak -t "${combchip}" -c "${combcontrol}" -f BED'" -n ${peakFile} -g ${GENOMESIZE} -p 1e-2 --nomodel --shiftsize ${fraglen} -B --SPMR" >> ${scriptName}
  fi
    echo 'rm -rf ${tmpdir}' >> ${scriptName}
    echo "rm -f ${peakFile}_peaks.xls ${peakFile}_peaks.bed ${peakFile}_summits.bed" >> ${scriptName}  
    # foldchange bedgraph
    if [[ ! -e ${fcFile} ]]
    then
	echo "macs2 bdgcmp -t ${peakFile}_treat_pileup.bdg -c ${peakFile}_control_lambda.bdg -o ${peakFile}.fc.bedgraph -m FE" >> ${scriptName}
	echo "slopBed -i ${peakFile}.fc.bedgraph -g ${CHRSIZE} -b 0 > ${peakFile}.fc.signal.bedgraph" >> ${scriptName}
	echo "bedGraphToBigWig ${peakFile}.fc.signal.bedgraph ${CHRSIZE} ${peakFile}.fc.signal.bw" >> ${scriptName}
	echo "gzip ${peakFile}.fc.signal.bedgraph" >> ${scriptName}
	echo "rm -f ${peakFile}.fc.bedgraph" >> ${scriptName}
    fi
    # LLR bedgraph
    if [[ ! -e ${llrFile} ]]; then
	echo "macs2 bdgcmp -t ${peakFile}_treat_pileup.bdg -c ${peakFile}_control_lambda.bdg -o ${peakFile}.llr.bedgraph -p 0.0001 -m logLR" >> ${scriptName}
	echo "slopBed -i ${peakFile}.llr.bedgraph -g ${CHRSIZE} -b 0 > ${peakFile}.llr.signal.bedgraph" >> ${scriptName}
	echo "bedGraphToBigWig ${peakFile}.llr.signal.bedgraph ${CHRSIZE} ${peakFile}.llr.signal.bw" >> ${scriptName}
	echo "gzip ${peakFile}.llr.signal.bedgraph" >> ${scriptName}
	echo "rm -rf ${peakFile}.llr.bedgraph" >> ${scriptName}
    fi
    echo "rm -f ${peakFile}_treat_pileup.bdg ${peakFile}_control_lambda.bdg" >> ${scriptName}
    
    # -------------------------
    # Submit script
    # -------------------------
    chmod 755 ${scriptName}
    cp ${scriptName} ${outFile}
    echo '======================================================================' >> ${outFile}
    if [[ ${MEM} -eq 0 ]]; then
	bsub -q research-rh6 -g "${JOBGRPID}" -J "${chipstub}" -W 24:00 -oo ${outFile} -eo ${errFile} < ${scriptName}
    else
	bsub -q research-rh6 -g "${JOBGRPID}" -J "${chipstub}" -W 24:00 -M "${memLimit}" -R "rusage[mem=${memLimit}]" -oo ${outFile} -eo ${errFile} < ${scriptName}
    fi
  
    # -------------------------
    # Delete temporary script
    # -------------------------	
    rm "${scriptName}"
    sleep 1s
done < "${IFILE:-/proc/${$}/fd/0}"