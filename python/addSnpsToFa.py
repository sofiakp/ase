import sys
import fileinput
import argparse
from random import randint
import os.path

def addSnps(fa, snppos, snpallele, p):
    for idx, pos in enumerate(snppos):
        alleles = snpallele[idx]
        if fa[pos].islower():
            fa[pos] = alleles[p].lower()
        else:
            fa[pos] = alleles[p]

def addSnpsX(mfa, pfa, snppos, snpallele):
    for idx, pos in enumerate(snppos):
        alleles = snpallele[idx]
        if mfa[pos].islower():        
            alleles = [a.lower() for a in alleles]
        if len(alleles) == 1:
            mfa[pos] = alleles[0]
        else:        
            mfa[pos] = alleles[1]
            pfa[pos] = alleles[0]

desc = 'Adds phased SNPs from a VCF to a genome to create two phased haplotypes'
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('fadir', help = 'Directory with reference fasta files. There should be one file chr*.fa for every chromosome.')
parser.add_argument('dict', help = 'List of chromosome names, one name per line. Output fasta files will have these chromosomes in the order specified in the file')
parser.add_argument('outpref', help = 'Output fasta files will be <outpref>.[pm]aternal.fa')
parser.add_argument('indiv', help = 'Sample name from the VCF for which haplotypes will be constructed')
parser.add_argument('-f', '--isfemale', help = 'If set to true, then chrY will be masked [%(default)s].', action = 'store_true', default = False)
parser.add_argument('-v', '--vcf', help = 'Vcf file. If not provided, reads from stdin', default = '')
parser.add_argument('-s', '--step', help = 'Number of characters per line of output fa [%(default)s].', type = int, default = 50)
parser.add_argument('-p', '--unphased', help = 'File to store unphased locations in VCF. These will be phased randomly.')
parser.add_argument('-c', '--chrom', help = 'Add chr to the chromosome names in VCF [%(default)s].', action = 'store_true', default = False)

args = parser.parse_args()
fadir = args.fadir
chromfile = args.dict
outpref = args.outpref
indiv = args.indiv
step = args.step
isfemale = args.isfemale
unphased = args.unphased
if len(unphased) > 0:
    unfile = open(unphased, 'w')

if args.vcf == '':
    vcflist = []
else:
    vcflist = [args.vcf]

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
        if args.chrom:
            chrom = "chr" + chrom
            if chrom == "chrMT":
                chrom = "chrM"
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
            print >> sys.stderr, 'Multiallelic or non-SNP, skipping', chrom, pos + 1
            continue
        subfields = fields[individx].split(':')
        genstr = subfields[0].split('|')
        if len(genstr) != 2:
            tmp = subfields[0].split('/')  
            if len(tmp) == 2: # randomly phase the unphased ones
                if len(unphased) > 0:
                    unfile.write(chrom + '\t' + str(pos + 1) + '\n')
                randidx = randint(0, 1)
                genstr = [tmp[randidx], tmp[1 - randidx]]
            elif len(tmp) == 1: # chrY
                genstr = tmp[0]
            else:
                print >> sys.stderr, 'Unknown genotype, skipping', chrom, pos + 1
                continue
        if genstr[0] not in ['0', '1'] or (len(genstr) > 1 and  genstr[1] not in ['0', '1']):
            print >> sys.stderr, 'Unknown genotype', chrom, pos + 1
            continue
        gen = [int(g) for g in genstr]
        snppos[chrom].append(pos)
        snpallele[chrom].append([alleles[g] for g in gen])

if len(unphased) > 0:
    unfile.close()

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
        if isfemale:
            pfa = ['N' for i in range(len(pfa))]
        else:
            addSnps(pfa, snppos[c], snpallele[c], 0)
    #elif c == 'chrX' and not isfemale:
    #    pfa = ['N' for i in range(len(pfa))]
    #    addSnps(mfa, snppos[c], snpallele[c], 1)
    elif c == 'chrX':
        addSnpsX(mfa, pfa, snppos[c], snpallele[c])
    else:
        addSnps(pfa, snppos[c], snpallele[c], 0)
        addSnps(mfa, snppos[c], snpallele[c], 1)


    for i in range(0, len(mfa), step):
        start = i
        stop = min(i + step, len(mfa))
        mout.write(''.join(mfa[start:stop]) + '\n')
        pout.write(''.join(pfa[start:stop]) + '\n')
mout.close()
pout.close()
