#!/bin/bash

PARDIR=$1
echo $PARDIR
while read -r indiv; do
    echo $indiv
    if [[ $indiv =~ ^NA ]]; then
	indiv=${indiv/NA/GM}
    fi
    indels=${PARDIR}/indels/${indiv}.indel.bed
    dels=${PARDIR}/sv/${indiv}.deletions.bed
    if [[ $indiv =~ ^GM12 ]]; then
	dups=${PARDIR}/sv/CEU.trio.tandem.bed
    else
	dups=${PARDIR}/sv/YRI.trio.tandem.bed
    fi
    echo "Merging $indels, $dels, $dups"
    cat $indels $dels $dups | sort -V | mergeBed -i stdin | sort -V > ${PARDIR}/${indiv}.variants.bed
done