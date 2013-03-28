import sys
import gzip
import fileinput
import re
import argparse

def compareChroms(fa1, fa2, chr1, chr2, snppos):
    print >> sys.stderr, 'Comparing', chr1, 'and', chr2
    if chr1 != chr2:
        print >> sys.stderr, 'Chromosomes were not in the same order'
        return
    if not chr1 in snppos:
        print >> sys.stderr, 'Skipping', chr1
        return
    if len(fa1) != len(fa2):
        print >> sys.stderr, 'Diff lengths for', chr1
        return
    
    for p in range(len(fa1)):
        if fa1[p] != fa2[p]:
            if p in snppos[chr1]:
                print chr1, p + 1, snppos[chr1][p]
            else:
                print chr1, p + 1, 'unexplained'

desc = 'Prints out SNP differences between two genomes. Assumes the chr lengths are the same.'
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('fa1', help = 'First genome')
parser.add_argument('fa2', help = 'Second genome. Chromosomes must be in the same order in the two genomes.')
parser.add_argument('dict', help = 'List of chromosome names, one name per line. SNPs not in the listed chromosomes will be ignored.')
parser.add_argument('indiv', help = 'Sample name from the VCF')
parser.add_argument('-v', '--vcf', help = 'Vcf file. If not provided, reads from stdin', default = '')

args = parser.parse_args()
fafile1 = args.fa1
fafile2 = args.fa2
chromfile = args.dict
indiv = args.indiv

if args.vcf == '':
    vcflist = []
else:
    vcflist = [args.vcf]

# Read genome dictionary
chroms = []
snppos = {}
with open(chromfile, 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        chroms.append(fields[0])
        snppos[fields[0]] = {}

# Read VCF
individx = -1
warnchroms = []
for line in fileinput.input(vcflist):
    if line[0:2] == '##':
        continue
    fields = line.strip().split()
    if fields[0] == '#CHROM':
        for i, f in enumerate(fields):
            if f == indiv:
                individx = i
                break
        if individx < 0:
            print >> sys.stderr, 'Individual', indiv, 'not found in vcf'
            exit(1)
    else:
        chrom = fields[0]
        if not chrom in chroms:
            if not chrom in warnchroms:
                print >> sys.stderr, 'Chromosome', chrom, 'in VCF is not in dictionary, skipping'
                warnchroms.append(chrom) # only warn once per chromosome
            continue
        if len(fields) < individx + 1:
            print >> sys.stderr, 'VCF line too short, skipping', chrom, pos + 1
            continue
        pos = int(fields[1]) - 1 # Convert 1-based VCF pos to 0-based index       
        alleles = fields[3:5]
        if len(alleles[1]) != 1:
            snppos[chrom][pos] = 'multi'
            continue
        subfields = fields[individx].split(':')
        genstr = subfields[0].split('|')
        if len(genstr) != 2:
            genstr = subfields[0].split('/')
            if len(genstr) == 1:
                snppos[chrom][pos] = 'haploid'
            elif len(genstr) == 2:
                snppos[chrom][pos] = 'unphased'
            else:
                nppos[chrom][pos] = 'missing'
fa1 = open(fafile1, 'r')
chr1 = re.sub('>|_maternal|_paternal', '', fa1.readline().strip())
fa2 = open(fafile2, 'r')
chr2 = re.sub('>|_maternal|_paternal', '', fa2.readline().strip())
fa1str = []
fa2str = []

for line1 in fa1:
    if line1[0] == '>': # Reached a new chromosome in fa1
        # Read a chromosome from file2
        for line2 in fa2:
            if line2[0] == '>': 
                # Reached a new chromosome in fa2, compare the previous chromosomes
                compareChroms(fa1str, fa2str, chr1, chr2, snppos)
                chr2 = re.sub('>|_maternal|_paternal', '', line2.strip())
                fa2str = []
                break
            else:
                fa2str.extend(line2.strip())
                if len(fa2str) % 10000000 == 0:
                    print >> sys.stderr, 'Read', len(fa2str), 'lines from', chr2
        chr1 = re.sub('>|_maternal|_paternal', '', line1.strip())
        fa1str = []
    else:
        fa1str.extend(line1.strip())
        if len(fa1str) % 10000000 == 0:
            print >> sys.stderr, 'Read', len(fa1str), 'lines from', chr1

for line2 in fa2:
    fa2str.extend(line2.strip())
compareChroms(fa1str, fa2str, chr1, chr2, snppos)

fa1.close()
fa2.close()
