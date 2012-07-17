import sys
import fileinput

for idx, line in enumerate(fileinput.input()):
    if (idx + 1) %  4== 0:
        for c in line.strip():
            if ord(c) > 74:
                print 1
                exit(0)
            if ord(c) < 64:
                print 0
                exit(0)
raise Exception('Ambiguous format')
