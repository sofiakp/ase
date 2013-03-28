import sys
import fileinput
import argparse

desc = "Converts VCF from b36 to hg19 coordinates.\n\
Updates the reference in the header but does not remove redundant fields (such as GP (b37 positions))\
SNPs with no mapping to hg19 or which map to a different chr are removed.\
Does NOT check whether the reference allele is indeed reference!!!!!\
SNP (dbSNP) IDs are removed (because they correspond to an old dbSNP version)."

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--bed', default = '',
                    help = 'BED file with mapped coordinates. Otherwise mapped coordinates are read from GP field in input VCF')
args = parser.parse_args()
liftchr = []
liftpos = []
if len(args.bed) > 0:
    with open(args.bed, 'r') as infile:
        for line in infile:
            fields = line.strip().split()
            liftchr.append(fields[0])
            liftpos.append(int(fields[1]) + 1)   

readpos = 0
for line in fileinput.input([]):
    if line.startswith('##reference'):
        print '##reference=hg19'
    elif line.startswith('#'):
        print line.strip()
    else:
        fields = line.strip().split()
        origChr = fields[0]
        found = False # found mapping in hg19?
        if len(args.bed) > 0:
            chr = liftchr[readpos]
            pos = str(liftpos[readpos])
            readpos = readpos + 1
            found = True
        else:
            subfields = fields[7].split(';')
            for f in subfields:
                if f.startswith('GP='):
                    found = True
                    chr = (f.split('=')[1]).split(':')[0]
                    pos = (f.split('=')[1]).split(':')[1]
                elif f.startswith('MP'):
                    found = False
                    break
        if not found or chr != origChr:
            continue
        fields[0] = 'chr' + chr
        if fields[0] == 'chrMT':
            fields[0] = 'chrM'
        fields[1] = pos
        fields[2] = '.'
        print '\t'.join(fields)
