import sys
#import fileinput
#import argparse
import numpy
import gzip

#geno = 
#with open(args.sites, 'r') as infile:
#    for line in infile:
#        fields = line.strip().split()
#        geno[':'.join(fields[0:2])] = fields[5]

filename = sys.argv[1]
if filename.endswith('.gz'):
    infile = gzip.open(filename, 'r')
else:
    infile = open(filename, 'r')

for line in infile:
    if line[0] == '#':
        continue
    fields = line.strip().split()
    pos = ':'.join(fields[0:2])
    #pat_fields = numpy.array([int(s) for s in fields[4].split(',')])
    #pat_fields = [str(s) for s in list(pat_fields[0:3:1] + pat_fields[3::1])]
    #mat_fields = numpy.array([int(s) for s in fields[5].split(',')])
    #mat_fields = [str(s) for s in list(mat_fields[0:3:1] + mat_fields[3::1])]
    #amb_fields = numpy.array([int(s) for s in fields[6].split(',')])
    #amb_fields = [str(s) for s in list(amb_fields[0:3:1] + amb_fields[3::1])]
    pat_fields = fields[4].split(',')
    mat_fields = fields[5].split(',')
    amb_fields = fields[6].split(',')
    print '\t'.join(['\t'.join(pat_fields), '\t'.join(mat_fields), '\t'.join(amb_fields)])
    #if any(pat_fields + mat_fields + amb_fields > 0):
    #    print line.strip()

infile.close()
