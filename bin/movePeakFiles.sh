while read inputline ; do
    #echo $inputline

    control=$(echo ${inputline} | awk '{print $2}' | sed -r 's/\.bam/.bed.gz/g') # extract second column as control file name, replace .bam with .bed.gz
    controlstub=$(echo ${control} | sed -r -e 's/;.*$//g' -e 's/\.bed\.gz//g') # Use first file name if multiple separated by ; and remove the .bed.gz extension to generate controlstub
    # extract first column as ChIP file name, replace .bam with .bed.gz and sort filenames
    chipstub=$(echo ${inputline} | awk '{print $3}')
    chip=$(echo ${inputline} | awk '{print $1}' | sed -r 's/\.bam/.bed.gz/g' | sed 's/;/\n/g')
    for cv in `echo $chip`; do
    if echo ${chip} | grep -q ';'; then

    else
	chipstub=$(echo ${chip} | sed -r 's/\.bed\.gz//g') # If one chip file then use remove .bed.gz from the file name and use that as ChIPstub
    fi
    peakPref="${chipstub}_VS_${controlstub}"
    outFile="${ODIR}/logFiles/${chipstub}_VS_${controlstub}.out"
    errFile="${ODIR}/logFiles/${chipstub}_VS_${controlstub}.err"
    peakFile=
    fcFile="${peakFile}.fc.signal.bedgraph.gz"
    llrFile="${peakFile}.llr.signal.bedgraph.gz"

done 