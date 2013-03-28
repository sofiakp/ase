import sys
import os.path
import fileinput

def switchgeno(g):
    if g == '0':
        return '1'
    if g == '1':
        return '0'
    return g

# Validates that the REF/ALT alleles of the VCF agree with the given fasta file.
# The VCF must be sorted.
# Writes to stdout, can read from file or stdin

fapref = sys.argv[1]

prevchr = ''
genoidx = 9
for line in fileinput.input(sys.argv[2:]):
    if line[0] == '#':
        print line.strip()
        continue
    fields = line.strip().split()
    chr = fields[0]
    pos = int(fields[1])
    ref = fields[3]
    alt = fields[4]
    if len(alt) != 1:
        print >> sys.stderr, 'Multiallelic or non-SNP, skipping', chr, pos
        continue
    if chr != prevchr:
        fa = []
        with open(os.path.join(fapref, chr + '.fa'), 'r') as fafile:
            for faline in fafile:
                if faline[0] == '>':
                    continue
                fa.extend(faline.strip().upper())
    if fa[pos - 1] != ref.upper() and fa[pos - 1] != alt.upper(): # vcf is 1-based
        print >> sys.stderr, 'Inconsistent pos:', chr, pos
    elif fa[pos - 1] == alt.upper():
        # Switch ref and alt
        newfields = []
        for f in fields[genoidx:]:
            print >> sys.stderr, f
            subfields = f.split(':')
            splitsymb = '/'
            geno = subfields[0].split(splitsymb)
            if len(geno) < 2:
                splitsymb = '|'
                geno = subfields[0].split(splitsymb)
            if len(geno) == 2:
                newfields.append(':'.join([switchgeno(geno[0]) + splitsymb + switchgeno(geno[1]), 
                                           ':'.join(subfields[1:])]))
            elif len(geno) == 1:
                newfields.append(':'.join([switchgeno(geno[0]), ':'.join(subfields[1:])]))
            else:
                print >> sys.stderr, 'Uknonwn or missing genotype', chr, pos
                print line.strip()
        fields[4] = ref
        fields[3] = alt
        print >> sys.stderr, 'Switched alleles at pos:', chr, pos
        print '\t'.join(fields[0:genoidx]) + '\t' + '\t'.join(newfields)
    else:
        print line.strip()
    prevchr = chr
