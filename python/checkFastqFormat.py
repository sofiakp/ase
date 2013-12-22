import sys
import fileinput
import argparse

desc = 'Checks FASTQ format. Reads a FASTQ from STDIN. Outputs 0 if the scores are Phred+33 (Illumina 1.8+) and 1 if the scores are Phred+64 (Illumina 1.5+). Raises an exception if the format cannot be determined (since the two ranges overlap).'
parser = argparse.ArgumentParser(description = desc)
args = parser.parse_args()

for idx, line in enumerate(fileinput.input()):
    if (idx + 1) %  4 == 0:
        for c in line.strip():
            if ord(c) > 74:
                print 1
                exit(0)
            if ord(c) < 64:
                print 0
                exit(0)
raise Exception('Ambiguous format')
