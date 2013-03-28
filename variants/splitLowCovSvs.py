import sys
import fileinput
import argparse
import re
import os.path

parser = argparse.ArgumentParser(description = 'Splits a VCF into BED files for each individual.')
parser.add_argument('vcf', help = 'Input VCF file')
parser.add_argument('-o', '--outdir', help = 'Directory for output files [Dir of input file]', default = '')
parser.add_argument('-s', '--outsuf', help = 'Suffix of output files [%(default)s]', default = 'sites.bed')
parser.add_argument('-v', '--outvcf', help = 'Write VCF output [%(default)s]', default = False, action = 'store_true')
parser.add_argument('-a', '--all', help = 'Include all variants (even those not present in the individual) [%(default)s]', default = False, action = 'store_true')

args = parser.parse_args()
outdir = args.outdir
outsuf = args.outsuf
outvcf = args.outvcf

if outdir == '':
    outdir = os.path.dirname(args.vcf)

outfiles = []
startcol = 9
titlePrint = False
for line in fileinput.input([args.vcf]):
    if line[0:2] == '##':
        continue
    fields = line.strip().split()
    if fields[0] == '#CHROM':
        for idx, f in enumerate(fields[startcol:]):
            indiv = re.sub('^NA', 'GM', f)
            if not titlePrint:
                outfiles.append(open(os.path.join(outdir, ''.join([indiv, outsuf])), 'w'))
                titlePrint = True
                if outvcf:
                    outfiles[idx].write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + indiv + '\n')
        continue

    chrom = fields[0]
    if re.match('^[0-9]+$', chrom) or chrom == 'X' or chrom == 'Y':
        chrom = 'chr' + chrom
    elif chrom == 'MT':
        chrom = 'chrM'

    start = int(fields[1])
    if not outvcf:
        start = start - 1 # VCF to BED
    name = fields[2]  
    ref = fields[3]
    alt = fields[4]
    vartype = ''
    # Indels and SVs have their starts and ends specified differently
    if re.search('VT=SNP', fields[7]):
        end = start + 1
        vartype = 'VT=SNP'
    elif re.search('INDEL', fields[7]):
        if len(ref) > len(alt): # deletion
            start = start + 1 # The bp before the deletion is in the ALT
            end = start + len(ref) - len(alt)
        else: # insertion
            end = start + 1
        vartype = 'VT=INDEL'
    else:
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
            end = start + 1
        start = start + istart
        end = end + iend
        vartype = 'VT=SV'
    if name == '.':
        name = '_'.join([chrom, str(start), str(end)])
    for idx, f in enumerate(fields[startcol:]):
        genos = f.split(':')[0]
        geno = re.split('\||/', genos)
        if (args.all or any([g == '1' for g in geno])) and (fields[6] == 'PASS' or re.search('PASS', fields[7]) or re.search('PASS', f)):
            if outvcf:
                if len(geno) == 2:
                    # Do not output pseudoautosomal
                    outfiles[idx].write('\t'.join([chrom, str(start), name, ref, alt, fields[5], fields[6], vartype, 'GT', genos]) + '\n')
            else:
                outfiles[idx].write('\t'.join([chrom, str(start), str(end), name, '.', '+']) + '\n')

for file in outfiles:
    file.close()
