import sys
import argparse 
import fileinput
import re

parser = argparse.ArgumentParser(description = 'Add validation information to a vcf file from another vcf')
parser.add_argument('vcf2', help = 'File from where to get additional validation information')
args = parser.parse_args()

valid = {}
with open(args.vcf2, 'r') as infile:
    for line in infile:
        if line[0] == '#':
            continue
        fields = line.strip().split()
        names = fields[2].split('&')
        if re.search('VALIDATED', fields[7]):
            for n in names:
                valid[n] = True

for line in fileinput.input([]):
    if line[0] == '#':
        print line.strip()
    else:
        fields = line.strip().split()
        name = fields[2]
        if not re.search('VALIDATED', fields[7]) and name in valid:
            fields[7] = ';'.join([fields[7], 'VALIDATED'])
        print '\t'.join(fields)
