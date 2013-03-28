import sys
import fileinput
import argparse

parser = argparse.ArgumentParser(description = "Reads VCF with SVs from STDIN and converts to a BED")
args = parser.parse_args()

for line in fileinput.input():
    if line[0] == '#':
        continue
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    name = fields[2]
    istart = 0
    iend = 0
    end = -1
    subfields = fields[7].split(';')
    for f in subfields:
        parts = f.split('=')
        if parts[0] == 'CIPOS':
            istart = int(parts[1].split(',')[0])
        elif parts[0] == 'CIEND':
            iend = int(parts[1].split(',')[1])
        elif parts[0] == 'END':
            end = int(parts[1])
    if end < start:
        end = start
    start = str(start + istart - 1)
    end = str(end + iend)
    if name == '.':
        name = '_'.join([chrom, start, end])
    print '\t'.join([chrom, start, end, name, '.', '+'])
