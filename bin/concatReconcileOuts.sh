#!/bin/bash

INDIR=$1
echo -e "File\tTotal\tPaternal\tPaternal_prc\tMaternal\tMaternal_prc\tAmbiguous\tAmbiguous_prc"
for i in `ls ${INDIR}/*out`; do
    base=`basename $i`
    base=${base/.out/}
    tot=`egrep "* Total" $i | cut -f2`
    p=`egrep "Num in readgroup paternal" $i | cut -f2-3`
    m=`egrep "Num in readgroup maternal" $i | cut -f2-3`
    a=`egrep "Num in readgroup amb" $i | cut -f2-3`
    echo -e "${base}\t${tot}\t${p}\t${m}\t${a}"
done