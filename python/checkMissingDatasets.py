import sys
import fileinput
import os.path
import re

# Reads a list of files (eg. mapped bam files) and a metadata file about the datasets
# and prints out files that are NOT in the metadata.

files = {}
# Read list of files from STDIN
for line in fileinput.input([]):
    pref = re.sub('SNYDER_HG19_|_reconcile.dedup.bam', '', os.path.basename(line.strip()))
    files[pref] = False
print len(files)

# Read metadata file
with open(sys.argv[1], 'r') as infile:
    for line in infile:
        fields = line.strip().split()
        if fields[0] == 'NAME' or fields[8] == 'NA': 
            continue
        name = fields[0].upper()
        if not name.startswith('SNYDER'):
            name = 'GM' + name
        if name in files:
            files[name] = True
        
for f, ex in files.iteritems():
    if not ex:
        print f
