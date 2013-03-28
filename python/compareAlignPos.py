import sys
import re

def clipPos(pos, cigar):
    if re.match('[0-9]*S', cigar):
        clip = int(cigar.split('S')[0])
        pos = pos - clip
    return pos

def main():
    # Print all reads in in1 that are mapped in different starting
    # pos in in2 or are not in in2 at all.
    # Assumes the input files do not contain pairs or duplicate names.

    in1 = sys.argv[1]
    in2 = sys.argv[2] 
    
    pos2 = {}
    with open(in2, 'r') as infile:
        for line in infile:
            fields = line.strip().split()
            name = fields[0]
            pos2[name] = clipPos(int(fields[3]), fields[5])
            
    with open(in1, 'r') as infile:
        for line in infile:
            fields = line.strip().split()
            name = fields[0]
            pos = clipPos(int(fields[3]), fields[5])
            if name in pos2:
                if pos != pos2[name]:
                    print '\t'.join([name, str(pos), str(pos2[name])])
            else:
                print '\t'.join([name, str(pos), '*'])

if __name__ == '__main__':
    main()
