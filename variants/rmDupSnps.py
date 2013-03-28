import sys
import fileinput

# Reads a VCF file from STDIN and removes duplicate SNPs (eg. same pos different dbSNP ids).

last_pos = ''
last_line = ''
for line in fileinput.input():
    if line[0] == '#':
        print line.strip()
        continue
    fields = line.strip().split()
    pos = fields[0] + ':' + fields[1]
    if last_line != '':
        if pos != last_pos:
            print last_line
        else:
            # If the SNP position is the same as in the previous line, then make sure the reference allele is the same.
            last_fields = last_line.split('\t')
            if last_fields[3] != fields[3]:
                print >> sys.stderr, 'Different reference allele for', fields[0], fields[1]
    last_pos = pos
    last_line = line.strip()

print last_line
