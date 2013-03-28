import sys
import fileinput

# Gets a BED file with counts in the last column and sums up the counts
# for the same region (first 3 columns).

last = ''
for line in fileinput.input():
    fields = line.strip().split()
    region = '_'.join(fields[0:3])
    if last != region:
        if last != '':
            print '\t'.join(['\t'.join(last.split('_')), str(counts)])
        counts = 0
    if fields[-1] != '.':
        counts = counts + int(fields[-1])
    last = region
if last != '':
    print '\t'.join(['\t'.join(last.split('_')), str(counts)])
