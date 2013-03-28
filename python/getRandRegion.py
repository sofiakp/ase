import sys
import random
import fileinput

chr_lens = {}
with open(sys.argv[1], 'r') as infile:
    for line in infile:
        fields = line.strip().split()
        chr_lens[fields[0]] = int(fields[1])

for line in fileinput.input([]):
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    length = int(fields[2]) - start
    rand_start = random.randint(1, chr_lens[chrom] - length + 1)
    print '\t'.join([chrom, str(rand_start), str(rand_start + length)])
