import sys
import argparse
import re

parser = argparse.ArgumentParser(description = 'Gets a VCF with genotypes for an individual, and a VCF with sites and writes the genotypes and phasing for the individual in STDOUT. This is similar to splitLowCovSvs but the input is different.')
parser.add_argument('genofile', help = 'VCF file for one individual')
parser.add_argument('posfile', help = 'VCF file with sites (actually, only the columns chrom, pos, id, and ref are used)')
parser.add_argument('fasta', help = 'Maternal personal genome of the individual. This is required to make sure the unphased SNPs agree between the personal genome and the genotype in the output file')
args = parser.parse_args()

fafile = open(args.fasta, 'r')
fa_chr = fafile.readline().strip().strip('>')

phases = {}
last_chrom = ''
switched = 0
# Get the individual's genotype at each position. SNPs not in this file are assumed to be homozygous ref.
with open(args.genofile, 'r') as infile:
    for line in infile:
        if line[0] == '#':
            continue
        fields = line.strip().split()
        if fields[0] != last_chrom:
            # Read a whole chromosome from the fasta file
            print >> sys.stderr, 'Saw', fields[0], 'reading', fa_chr
            mfa = []
            for faline in fafile:
                if faline[0] == '>':
                    if not faline.startswith('>' + fields[0]):
                        fa_chr = faline.strip().strip('>')
                        break
                else:
                    mfa.extend(faline.strip().upper())
        last_chrom = fields[0]

        # If the genotype says that maternal is alt but maternal is actually ref, then switch the genotype
        if fields[9] == '0/1' and mfa[int(fields[1]) - 1] == fields[3]:
            fields[9] = '1/0'
            switched = switched + 1
            #print >> sys.stderr, 'Switching', fields[0], fields[1]
        elif fields[9] == '1/0' and mfa[int(fields[1]) - 1] != fields[3]:
            fields[9] = '0/1'
            switched = switched + 1
            #print >> sys.stderr, 'Switching', fields[0], fields[1]
        # Map the position of the SNP to the ref allele, alt allele, and genotype.
        phases[fields[0] + ':' + fields[1]] = '_'.join([fields[3], fields[4], fields[9]])
print >> sys.stderr, 'Total switched', switched

print '\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'GENO'])
# Read the list of sites. The output will contain one line for each site in this file.
with open(args.posfile, 'r') as infile:
    for line in infile:
        if line[0] != '#':
            fields = line.strip().split()
            pos = fields[0] + ':' + fields[1]
            if pos in phases:
                pfields = phases[pos].split('_')
                #isphased = int(pfields[2].find('/') < 0)
                #parts = re.split('\||/', pfields[2])
                print '\t'.join([fields[0], fields[1], fields[2], pfields[0], pfields[1], pfields[2]])
            else:
                print '\t'.join([fields[0], fields[1], fields[2], fields[3], fields[3], '0|0']) # hom ref and (trivially) phased
