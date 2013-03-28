import sys
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('chrLen', help = 'chromosome lengths')
parser.add_argument('-w', '--win', type = int, help = 'window length')
args = parser.parse_args()
win = args.win

with open(args.chrLen, 'r') as infile:
    for line in infile:
        fields = line.strip().split()
        if line[0] == '#' or fields[0] == 'chrom':
            continue
        chrom = fields[0]
        chr_len = float(fields[1])
        nchunks = int(math.ceil(chr_len / win))
        for s in range(nchunks):
            print '\t'.join([chrom, str(s * win), str(min(chr_len, (s + 1) * win))])

#outpref = sys.argv[1]
#chrmin = 20000000
#chrmax = 102531392
#chrlen = chrmax - chrmin

#wins = [5000, 10000, 50000, 100000, 1000000, 50000000]
#for i in range(len(wins)):
#    win = wins[i]
#    step = 5000
#    nwin = int(math.ceil(chrlen / step))
#    #outfile = fopen(outpref + '_w' + str(win) + '.bed', 'w')
#    for s in range(0, nwin):
#        start = chrmin + s * step
#        end = min(chrmin + chrlen, start + win)
#        print '\t'.join(['chr15', str(start), str(end), str(win)])
