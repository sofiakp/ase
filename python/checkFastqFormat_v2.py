import sys
import fileinput

minq = 10000
maxq = 0
minchar = ''
maxchar = ''
for idx, line in enumerate(fileinput.input()):
    if (idx + 1) %  4== 0:
        if idx > 10000:
            break
        for c in line.strip():
            if ord(c) < minq:
                minq = ord(c)
                minchar = c
            if ord(c) > maxq:
                maxq = ord(c)
                maxchar = c
print minq, minchar, maxq, maxchar
