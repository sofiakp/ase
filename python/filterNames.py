import sys
import fileinput

# Selects reads with the given names from a SAM file

namefile = sys.argv[1]

names = {}
with open(namefile, 'r') as infile:
    for line in infile:
        names[line.strip()] = True

for line in fileinput.input(sys.argv[2:]):
    if line[0] == '@':
        print line.strip()
    fields = line.strip().split()
    if fields[0] in names:
        print line.strip()
