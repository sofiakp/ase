import sys
import argparse
import os.path

def addSnps(fa, snppos, alleles):
    for idx, pos in enumerate(snppos):
        if fa[pos].islower():
            fa[pos] = alleles[idx].lower()
        else:
            fa[pos] = alleles[idx]

def writeChrom(c, mout, pout, mfa, pfa, step):
    if c == 'chrY':
        mfa = ['N' for i in range(len(mfa))]
        if isfemale:
            pfa = ['N' for i in range(len(pfa))]

    mout.write('>' + c + '\n')
    pout.write('>' + c + '\n')
    for i in range(0, len(mfa), step):
        start = i
        stop = min(i + step, len(mfa))
        mout.write(''.join(mfa[start:stop]) + '\n')
        pout.write(''.join(pfa[start:stop]) + '\n')

desc = 'Adds phased SNPs from Seqphase to a genome to create two phased haplotypes'
parser = argparse.ArgumentParser(description = desc)
parser.add_argument('fadir', help = 'Directory with reference fasta files. There should be one file chr*.fa for every chromosome.')
parser.add_argument('dict', help = 'List of chromosome names, one name per line. Output fasta files will have these chromosomes in the order specified in the file')
parser.add_argument('mappref', help = 'Prefix of map and ped files. These should be named <map>_chrX.[map|ped].')
parser.add_argument('outpref', help = 'Output fasta files will be <outpref>.[pm]aternal.fa')
parser.add_argument('indiv', help = 'Sample name.')
parser.add_argument('-f', '--isfemale', help = 'If set to true, then chrY will be masked [%(default)s].', action = 'store_true', default = False)
parser.add_argument('-p', '--seqphase', help = 'Prefix of seqphase files. There should be one file <seqphase>chrX.txt for each chromosome.', 
                    default = '')
parser.add_argument('-s', '--step', help = 'Number of characters per line of output fa [%(default)s].', type = int, default = 50)
#parser.add_argument('-c', '--chrom', help = 'Add chr to the chromosome names in VCF [%(default)s].', action = 'store_true', default = False)

args = parser.parse_args()
fadir = args.fadir
outpref = args.outpref
indiv = args.indiv
step = args.step
isfemale = args.isfemale

# Read genome dictionary
chroms = []
with open(args.dict, 'r') as infile:
    for line in infile:
        fields = line.strip().split('\t')
        chroms.append(fields[0])

mout = open(outpref + '.maternal.fa', 'w')
pout = open(outpref + '.paternal.fa', 'w')
for c in chroms:
    mfa = []
    pfa = []
    snppos = []
    with open(os.path.join(fadir, c + '.fa'), 'r') as fafile:
        print >> sys.stderr, 'Loading', c
        for faline in fafile:
            if faline[0] != '>':
                mfa.extend(faline.strip())
                pfa.extend(faline.strip())    
                
    if os.path.isfile(args.mappref + c + '.ped'):
        # Read ped file to get the index of the individual
        individx = -1
        with open(args.mappref + c + '.ped', 'r') as infile:
            for idx, line in enumerate(infile):
                if line.strip().split()[0] == 'SNYDER_HG19_' + indiv:
                    individx = idx
        assert(individx > -1)
        
        with open(args.mappref + c + '.map', 'r') as infile:
            for line in infile:
                fields = line.strip().split()
                snppos.append(int(fields[3]) - 1) # Convert to 0-based Pythong index
        print >> sys.stderr, 'Read', str(len(snppos)), 'SNPS'

        with open(args.seqphase + c + '.txt', 'r') as genotfile:
            for idx, line in enumerate(genotfile):
                if idx == 2 * individx:
                    alleles = line.strip().split()
                    assert all([a in ['A', 'C', 'G', 'T'] for a in alleles])
                    assert len(snppos) == len(alleles)
                    addSnps(pfa, snppos, alleles)
                elif idx == 2 * individx + 1:
                    alleles = line.strip().split()
                    assert all([a in ['A', 'C', 'G', 'T'] for a in alleles])
                    assert len(snppos) == len(alleles)
                    addSnps(mfa, snppos, alleles)
                elif idx > 2 * individx + 1:
                    break
    else:
        print >> sys.stderr, 'Missing MAP/PED files for', c

    writeChrom(c, mout, pout, mfa, pfa, step)

mout.close()
pout.close()
