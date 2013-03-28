import sys
import fileinput

#sam = open(sys.argv[1], 'r')
outpref = sys.argv[1]
common = open(outpref + 'common.pairs.txt', 'w')
bad = open(outpref + 'bad.pairs.txt', 'w')
bestm = open(outpref + 'maternal.pairs.txt', 'w')
bestp = open(outpref + 'paternal.pairs.txt', 'w')

prevname = ''
firstmat = (-1, 'X')
secmat = (-1, 'X')
firstpat = (-1, 'X')
secpat = (-1, 'X')
for line in fileinput.input(sys.argv[2:]):
    fields = line.strip().split()
    name = fields[0]
    if len(prevname) > 0 and name != prevname:
        # Pair wasn't found in paternal file
        if firstpat[0] < 0:
            if secpat[0] >= 0:
                print >> sys.stderr, prevname, 'is unmated in father'
            else:
                bestm.write(prevname + '\n')
        elif firstmat[0] < 0:
            if secmat[0] >= 0:
                print >> sys.stderr, prevname, 'is unmated in mother'
            else:
                bestp.write(prevname + '\n')
        else:
            qualm = max(firstmat[0], secmat[0])
            qualp = max(firstpat[0], secpat[0])
            if qualm < qualp:
                bestp.write(prevname + '\n')
            elif qualm > qualp:
                bestm.write(prevname + '\n')
            elif firstmat[1] == firstpat[1] and secmat[1] == secpat[1]:
                # equal qualities, same positions
                common.write(prevname + '\n')
            else:
                # soft-clipping is ignored, so the actual mapping position
                # might still be the same
                bad.write(prevname + '\n')
        firstmat = (-1, 'X')
        secmat = (-1, 'X')
        firstpat = (-1, 'X')
        secpat = (-1, 'X')
    flag = bin(int(fields[1]))
    if len(flag) >  8: # The first two chars are 0b
        isfirst = bool(flag[-7] == '1')
    else:
        isfirst = False
    qual = int(fields[4]) #-10**(-int(fields[4])/10)
    chr = fields[2].split('_')[0]
    pos = fields[3]
    ispat = fields[2].find('_p') >= 0
    tuple = (qual, chr + ':' + pos)
    if ispat:
        if isfirst:
            firstpat = tuple
        else:
            secpat = tuple
    else:
        if isfirst:
            firstmat = tuple
        else:
            secmat = tuple
    prevname = name
    
#sam.close()
common.close()
bad.close()
bestm.close()
bestp.close()
