import sys
import fileinput
import argparse
from random import randint
import os.path
import re

def addSnps(fa, snppos, snpallele, p):
    for idx, pos in enumerate(snppos):
        alleles = snpallele[idx]
        if fa[pos].islower():
            fa[pos] = alleles[p].lower()
        else:
            fa[pos] = alleles[p]

desc = 'Adds Snyder SNPs to genome based on a custom-format file read from STDIN'
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('fadir', help = 'Directory with reference fasta files. There should be one file chr*.fa for every chromosome.')
parser.add_argument('dict', help = 'List of chromosome names, one name per line. Output fasta files will have these chromosomes in the order specified in the file')
parser.add_argument('outpref', help = 'Output fasta files will be <outpref>.[pm]aternal.fa')
parser.add_argument('-s', '--step', help = 'Number of characters per line of output fa [%(default)s].', type = int, default = 50)
args = parser.parse_args()
fadir = args.fadir
chromfile = args.dict
outpref = args.outpref
step = args.step

# Read genome dictionary
chroms = []
snppos = {}
snpallele = {}
with open(chromfile, 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        chroms.append(fields[0])
        snppos[fields[0]] = []
        snpallele[fields[0]] = []

# Read SNP file
warnchroms = []
for line in fileinput.input([]):
    fields = line.strip().split()
    if fields[0] == 'chr':
        continue
    chrom = fields[0]
    if not chrom in chroms:
        if not chrom in warnchroms:
            print >> sys.stderr, 'Chromosome', chrom, 'in VCF is not in dictionary, skipping'
            warnchroms.append(chrom) # only warn once per chromosome
        continue
    pos = int(fields[2]) - 1 # Convert 1-based VCF pos to 0-based index
    if re.match('^[ACGT]$', fields[5]) and re.match('^[ACGT]$', fields[6]):
        snppos[chrom].append(pos)
        snpallele[chrom].append(fields[5:7])

mout = open(outpref + '.maternal.fa', 'w')
pout = open(outpref + '.paternal.fa', 'w')
for c in chroms:
    mfa = []
    pfa = []
    with open(os.path.join(fadir, c + '.fa'), 'r') as fafile:
        print >> sys.stderr, 'Loading', c
        for faline in fafile:
            if faline[0] == '>':
                continue
            mfa.extend(faline.strip())
            pfa.extend(faline.strip())
    mout.write('>' + c + '\n')
    pout.write('>' + c + '\n')
 
    # Mask both maternal and paternal chrY for females
    # and maternal one for males.
    # Don't mask chrX for males because you will lose
    # the pseudoautosomal region.
    if c == 'chrY':
        mfa = ['N' for i in range(len(mfa))]
        addSnps(pfa, snppos[c], snpallele[c], 1)
    else:
        addSnps(pfa, snppos[c], snpallele[c], 1)
        addSnps(mfa, snppos[c], snpallele[c], 0)

    print len(mfa), len(pfa)
    for i in range(0, len(mfa), step):
        start = i
        stop = min(i + step, len(mfa))
        mout.write(''.join(mfa[start:stop]) + '\n')
        pout.write(''.join(pfa[start:stop]) + '\n')
mout.close()
pout.close()
