import sys
import fileinput 
import argparse
import re 
import os.path

parser = argparse.ArgumentParser(description = 'Splits a VCF with SVs into BED files for each individual.')
parser.add_argument('mapped', help = 'BED file with the mapped variants')
parser.add_argument('-o', '--outdir', help = 'Directory for output files [Dir of input file]', default = '')
parser.add_argument('-s', '--outsuf', help = 'Suffix of output files [%(default)s]', default = 'sites.bed')
parser.add_argument('-v', '--validOnly', help = 'Only select VALIDATED variants', action = 'store_true', default = False)
parser.add_argument('-t', '--sitesOnly', help = 'Do not split by individual, just report the sites', action = 'store_true', default = False)
parser.add_argument('-p', '--passOnly', help = 'Search for flag PASS for each individual and output these sites only', 
                    action = 'store_true', default = False)
args = parser.parse_args()
outdir = args.outdir
outsuf = args.outsuf
validOnly = args.validOnly
sitesOnly = args.sitesOnly
passOnly = args.passOnly

if outdir == '':
    outdir = os.path.dirname(args.mapped)

mapped = {}
with open(args.mapped, 'r') as infile:
    for line in infile:
        fields = line.strip().split()
        fields[0] = "chr" + fields[0]
        # Mark as unmapped those that have multiple mappings
        if int(fields[4]) > 1:
            mapped[fields[3]] = ''
        else:
            mapped[fields[3]] = fields[0:4]

outfiles = []
startcol = 9
for line in fileinput.input([]):
    if line[0:2] == '##':
        continue
    fields = line.strip().split()
    if fields[0] == '#CHROM':
        if sitesOnly:
            outfiles.append(open(os.path.join(outdir, outsuf), 'w'))
        else:
            for f in fields[startcol:]:
                indiv = re.sub('^NA', 'GM', f)
                outfiles.append(open(os.path.join(outdir, ''.join([indiv, outsuf])), 'w'))
    else:
        info = fields[7]
        if validOnly and not re.search('VALIDATED', info):
            continue
        name = fields[2]
        if name == '.':
            name = '_'.join([fields[0], str(int(fields[1]) - 1), fields[1]])
        if not name in mapped or mapped[name] == '':
            continue
        region = mapped[name]
        if sitesOnly:
            outfiles[0].write('\t'.join(region) + '\n')
        else:
            for idx, f in enumerate(fields[startcol:]):
                parts = f.split(':')
                geno = re.split('\||/', parts[0])
                if any([g == '1' for g in geno]) and (not passOnly or re.search('PASS', f)):
                    outfiles[idx].write('\t'.join(region) + '\n')

for file in outfiles:
    file.close()
