#!/bin/bash
# Takes one argument which is the file containing COLUMN1=tagAlign file names and COLUMN2 the corresponding extension lengths

# ###########################################
# AUXILIARY FUNCTIONS
# ###########################################

# Checks if a file or directory exists and exits if it fails
FUNCdoesFileExist()
{
    local VarName=$1
    local fileName=$2    
    if [[ ! -e "${fileName}" ]]
    then
        echo "ERROR in variable ${VarName}:${fileName} - File/Directory does not exist"
        exit 1
    fi
}

# ###########################################
# MAIN
# ###########################################

if [[ "$#" -lt 4 ]]
then
	echo "$(basename $0): Will create signal files from a set of tagAlign/BAM files and corresponding extension lengths" >&2
	echo "USAGE: $(basename $0) <PairFile> <FraglenFile> <inputDir> <outputDir> <replaceFlag> <memory> <oformat> <kernelName> <normFlag>" >&2
	echo ' <PairFile>: Col1: ChIP replicate file names separated by ; Col2: OutputFilePrefix' >&2
	echo ' <FraglenFile>: tab delimited file where Col1: path_to_file and Col2: fragment length and Col3: contains smoothing window size' >&2
	echo ' <inputDir>: Root directory under which all the align files are located' >&2
	echo ' <outputDir>: directory where you want the output files to be' >&2
	echo ' <replaceFlag>: TRUE/FALSE If set to TRUE then files will be overwritten else skipped' >&2
	echo ' <memory>: in GB' >&2
	echo ' <oformat>: mat,bg,wig' >&2
	echo ' <kernelName>: rectangular,triangular,epanechnikov,biweight,triweight,cosine,gaussian,tukey' >&2
	echo ' <normFlag>: 0: counts, 5: foldchange' >&2
	echo ' <mapFilter>: 0 (allows all locations)' >&2
	exit 1
fi

PAIRFILE=$1 # File containing ChIP file names and output file prefix
FUNCdoesFileExist 'PAIRFILE' "${PAIRFILE}"

FRAGLENFILE=$2 # File containing fragment lengths and smoothing sizes
FUNCdoesFileExist 'FRAGLENFILE' "${FRAGLENFILE}"

IDIR=$3 # Input root directory containing tagAlign/BAM files
FUNCdoesFileExist 'IDIR' "${IDIR}"
printf "ALIGNDIR: %s\n" ${IDIR}

ODIR=$4 # Output directory
FUNCdoesFileExist 'ODIR' "${ODIR}"
printf "OUTPUTDIR: %s\n" ${ODIR}

REPLACEFLAG=$5 # If set to TRUE then if output file exists, it is replaced. Else it is skipped
printf "REPLACEFLAG: %s\n\n" ${REPLACEFLAG}

MEM=$6

OFORMAT=$7

MAXJOBS=30

KERNELNAME=$8

JOBGROUP="/wiggler${RANDOM}"
MAXJOBS=30
#bgadd -L ${MAXJOBS} $JOBGROUP

NORMFLAG=$9 # Normalization flag
printf "NORMFLAG: %d\n" ${NORMFLAG}

MAPFILTER=${10}
printf "MAPFILTER: %s\n" ${MAPFILTER}

printf "OUTPUTFORMAT: %s\n" ${OFORMAT}
if [[ "${OFORMAT}" == 'bg' ]]
then
	OFEXT=".norm${NORMFLAG}.rawsignal.bedgraph" # Output file extension
else
	OFEXT=".norm${NORMFLAG}.rawsignal.${OFORMAT}" # Output file extension
fi
LOGEXT=".norm${NORMFLAG}.rawsignal.log" # log file extension

CHR_DIR="${SEQDIR}/encodeHg19Male"

UNIQ_DIR="${UMAPDIR}/hg19_chromvar_k101/globalmap_k1tok1000"

BASECMD="align2rawsignal -mm=${MEM} -of=${OFORMAT} -n=${NORMFLAG} -f=${MAPFILTER} -k=${KERNELNAME} -s=${CHR_DIR} -u=${UNIQ_DIR}" # Base ta2rs command                                                                                                                        

count=0                                                                                             
while read line 
do
    chipfiles=$(echo ${line} | awk '{print $1}')
    outputPrefix=$(echo ${line} | awk '{print $2}')

    # ====================================	    
    # Generate logfile name
    # ====================================
    LOGFILE="${ODIR}/${outputPrefix}${LOGEXT}"    
    OFNAME="${ODIR}/${outputPrefix}${OFEXT}"
    CMD="$BASECMD -v=${LOGFILE}" # update command with --v
    
    # ====================================   
    # Get files and tag extensions corresponding to the root
    # ====================================
    nf=0
    for i in $(echo ${chipfiles} | sed -r 's/;/\n/g')
    do
	c=$((c+1))
	filepath=$(find ${IDIR} -name $i | head -1)
	if [[ ! -e ${filepath} ]]
	then
	    echo "WARNING: File $i not found" >&2
	    nf=1
	    continue
	else
	    fragline=$(grep -F $i ${FRAGLENFILE} | head -1)
	    if [[ -z ${fragline} ]]
	    then
		echo "WARNING: No fragment length information for file $i" >&2
		nf=1
	    else
		fraglen=$(echo ${fragline} | awk '{print $2}')
		smoothlen=$(echo ${fragline} | awk '{print $3}')
		CMD="$CMD -i=${filepath} -l=${fraglen}"
	    fi
	fi
    done
    CMD="$CMD -w=${smoothlen}"
    if [[ ${nf} -eq 1 ]]
    then
        echo "Error for ${line}" >&2
	continue
    fi
    
    # ====================================
    # Generate outputfile name
    # ====================================
    # Update command with --ofile and optionally --output-max-tags
    if [[ "${OFORMAT}" == 'mat' ]]
    then
        CMD="$CMD -o=${OFNAME}" # Update command with --ofile
	#MAXTAGSFNAME=$(echo ${OFNAME} | sed -r 's/\.rawsignal\./.rawsignal.maxtags./g')
    	#CMD="$CMD --output-max-tags=${MAXTAGSFNAME}"
    else
    	OFNAME="${OFNAME}.gz"
    	CMD="$CMD | grep -v Warning | gzip -c 1> ${OFNAME}"
    fi

    # ====================================
    # Check if output file exists and if REPLACEFLAG is not set then skip
    # ====================================
    if [[ -e "${OFNAME}" && "${REPLACEFLAG}" != 'TRUE' ]]
    then
        echo "Output File Exists: Skipping ${OFNAME}"
        printf "%s\t%s\tSkipped\n" $(date +%D) ${OFNAME} >> skippedDatasets.log
        continue
    fi

    # Create temp submit/run script
    SUBMITFILE="${outputPrefix}.sh"
    echo '#!/bin/bash' > "$SUBMITFILE"
    echo "TEMPDIR=${TMP}/"'temp_${RANDOM}${RANDOM}' >> "$SUBMITFILE" # Create temp directory for source files
    echo 'mkdir $TEMPDIR' >> "$SUBMITFILE"
    echo 'export MCR_CACHE_ROOT=${TEMPDIR}' >> "$SUBMITFILE"
    echo "cp $(which align2rawsignal) "'${TEMPDIR}/' >> "$SUBMITFILE"
    echo 'sleep 10s' >> "$SUBMITFILE"
    echo 'cd $TEMPDIR' >> "$SUBMITFILE"
    echo 'sleep 2s' >> "$SUBMITFILE"
    #echo "export TMP=${ENSEMBL_SCRATCH}/tmp" >> "$SUBMITFILE" # export TMP to ENSEMBL_TMP
    echo "./${CMD}" >> "$SUBMITFILE"
    echo 'cd ..' >> "$SUBMITFILE"
    echo 'sleep 2s' >> "$SUBMITFILE"
    echo 'rm -rf $TEMPDIR' >> "$SUBMITFILE"
    chmod 755 "${SUBMITFILE}"

    count=$(( count + 1 ))
    if [[ $(( count % MAXJOBS )) -eq 0 ]]
	then
	echo 'Waiting ..'
	sleep 10m
    fi

    ALLOC=$(( MEM * 1024 ))
    # Submit/Run the script
    bsub -J ${outputPrefix} -g "${JOBGROUP}" -M ${ALLOC} -R "rusage[mem=${ALLOC}]" -o "${LOGFILE}" -e "${LOGFILE}.err" < "${SUBMITFILE}"
    sleep 2s

    rm "${SUBMITFILE}"
done < ${PAIRFILE}
