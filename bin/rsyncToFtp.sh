#!/bin/bash

opts="-rtlgpuav"
ftpRoot=${FTPDIR}/temp/chromatinVariation
data=${MAYAROOT}/rawdata

######### ALWAYS TEST WITH -n AFTER MAKING CHANGES
rsync $opts --delete --exclude="*.RData" --exclude="rdata" --exclude="fc/matFiles/" --exclude="llr/matfiles" ${data}/signal $ftpRoot
rsync $opts --delete --exclude="*.RData" --exclude="rdata" --exclude="matFiles" ${data}/genomeGrid $ftpRoot
rsync $opts --delete --exclude="*.RData" --exclude="rdata" --exclude="matFiles" --exclude="gencodev13.noM" ${data}/transcriptomes $ftpRoot
rsync $opts ${data}/metadata $ftpRoot