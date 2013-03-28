import sys
import fileinput
import argparse
import re

parser = argparse.ArgumentParser(description = 'Reads a VCF from STDIN and selects variants that are heterozygous for the given sample.')
parser.add_argument('-i', '--indiv', help = 'Individual to select. Default: first individual', default = '')
args = parser.parse_args()
sample = args.indiv

startcol = 9
sampleIdx = -1
for line in fileinput.input([]):
    if line.startswith('##'):
        print line.strip()
        continue
    fields = line.strip().split()
    if fields[0] == '#CHROM':
        for idx, f in enumerate(fields[startcol:]):
            if f == sample:
                sampleIdx = idx + startcol
                break
            if sampleIdx < 0:
                sampleIdx = startcol
        print '\t'.join(fields[0:startcol]) + '\t' + fields[sampleIdx]
    else:
        geno = re.split('\||/', fields[sampleIdx])
        if len(geno) == 2 and geno[0] != geno[1] and all([g in ['0', '1'] for g in geno]):
            print '\t'.join(fields[0:startcol]) + '\t' + fields[sampleIdx]
