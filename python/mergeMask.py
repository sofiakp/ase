import sys

def outputBEDline(chrom, startCord, endCord):
    outs = 'chr' + str(chrom) + '\t' + str(startCord) + '\t' + str(endCord) + '\n'
    outfile.write(outs)

outfile = open('mergedMask.bed','w')
chromNum = range(1, 23)
chromNum.append('X')
effMask = {'L', 'H', 'Z', 'Q'}
for chrom in chromNum:
    print '* On chromesome {}.'.format(chrom)
    fileName = 'chr' + str(chrom) + '.pilot_style_mask.fasta'
    infile = open(fileName)
    infile.readline()
    startCord = -1
    endCord = -1
    cord = -1
    formerChar = '#'
    for line in infile:
        for j in range(len(line) - 1):
            cord += 1
            if formerChar == '#':
                if line[j] in effMask:
                    formerChar = line[j]
                    startCord = cord
                else:
                    continue
            else:
                if line[j] == formerChar:
                    continue
                else:
                    if line[j] in effMask:
                        continue
                    else:
                        formerChar = '#'
                        endCord = cord
                        outputBEDline(chrom, startCord, endCord)
    if formerChar != '#':
        outputBEDline(chrom, startCord, cord + 1)
