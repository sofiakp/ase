import sys
import fileinput
import os.path
import re

def selectField(f1, f2):
    if len(f1) == 0 and len(f2) == 0:
        return ''
    elif len(f1) > len(f2):
        return f1
    elif len(f1) < len(f2):
        return f2
    else:
        if f1 == f2:
            return f1
        else:
            return ''

# Create a text file with sample info from Maya's google docs
# Docs are (assumed to be) downloaded as CSV
# Reads from file or stdin

print '\t'.join(['Individual', 'Mark', 'Replicate', 'Sample_name',
                 'Read1', 'Read2', 'http_location'])

names = {}
for line in fileinput.input():
    if line.find('http') < 0: # heuristic for removing all header lines
        continue
    else:
        fields = [f.strip() for f in line.strip().split(',')]
        nameFields = [f.strip() for f in fields[0].split('_')] 
        # Mark is included in two different fields. Pick "best".
        # Remove '.' (e.g. H2A.Z is sometimes written as H2AZ)
        mark = selectField(fields[4].replace('.', ''), 
                           nameFields[1].replace('.', ''))
        if len(mark) == 0:
            print >> sys.stderr, 'Mark missing or ambiguous', line.strip()
            exit(1)
        # Make some effort to have consistent names
        # mark = re.sub('ac|AC|aC', 'Ac', mark)
        # mark = re.sub('ME|Me|mE', 'me', mark)
        # if mark.upper() == 'INPUT':
        #    mark = 'Input'
        mark = mark.upper()

        indiv = selectField(fields[3], nameFields[0])
        if len(indiv) == 0:
            print >> sys.stderr, 'Individual name missing or ambiguous', line.strip()
        # Sometimes mark, indiv, and replicate don't uniquely identify sample
        addfields = ''
        if len(nameFields) < 3:
            rep = '1' # replicate
        else:
            rep = nameFields[2]
            if not rep.isdigit():
                rep = '1'
                addfields = '_'.join(nameFields[2:])
            elif len(nameFields) > 3:
                # TR at the end of name...
                addfields = '_'.join(nameFields[3:])

        sampleName = '_'.join([mark, indiv, rep, addfields]).strip('_')

        # In the future add some extra check to ignore fq files that appeared in
        # earlier sample lists (from previous downloads) and check that 
        # new names generated have not appeared before.
        if sampleName in names:
            print >> sys.stderr, 'Duplicate sample name', sampleName
            initName = sampleName
            while sampleName in names:
                names[initName] = names[initName] + 1
                sampleName = initName + str(names[initName]) 
        else:
            names[sampleName] = 1
        read1loc = fields[8]
        read2loc = fields[9]
        (loc1, fqpref1) = os.path.split(read1loc)
        (loc2, fqpref2) = os.path.split(read2loc)
        if loc1 != loc2:
            print >> sys.stderr, 'Different http locations for the two fq files', \
                loc1,loc2
            exit(1)

        outfields = [indiv, mark, rep, sampleName, fqpref1, fqpref2, loc1]
        if all([len(f) > 0 for f in outfields]):
            print '\t'.join(outfields)
        else:
            print >> sys.stderr, 'Info missing for line', line.strip()
