import sys
import os.path

def outputBEDline(chrom, startCord, endCord):
    outs = 'chr' + str(chrom) + '\t' + str(startCord) + '\t' + str(endCord) + '\n'
    print outs

indir = sys.argv[1]
chromNum = range(1, 23)
chromNum.append('X')
effMask = ['N'] #{'L', 'H', 'Z', 'Q'}
for chrom in chromNum:
    print >> sys.stderr, '* On chromosome', chrom
    fileName = os.path.join(indir, 'chr' + str(chrom) + '.fa')
    infile = open(fileName)
    infile.readline()
    startCord = -1
    endCord = -1
    cord = -1
    formerChar = '#'
    for line in infile:
        for j in range(len(line) - 1):
            cord += 1
            if line[j] in effMask:
                #print 'chr' + str(chrom) + '\t' + str(cord)
                if formerChar == '#':
                    #if line[j] in effMask:
                    formerChar = line[j]
                    startCord = cord
            else:
                #if line[j] == formerChar:
                #    continue
                #else:
                #    if line[j] in effMask:
                #        continue
                #    else:
                formerChar = '#'
                endCord = cord
                if startCord > -1:
                    outputBEDline(chrom, startCord, endCord)
    if formerChar != '#':
        outputBEDline(chrom, startCord, cord + 1)
