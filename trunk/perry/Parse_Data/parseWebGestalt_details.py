#!/usr/bin/python
import sys

if len(sys.argv) <= 3:
    print 'enter the fileName & the order you want the fields in'
    sys.exit(0)

def removeQuotes(s):
    return s[1:-2].strip()

order = []
for o in sys.argv[2:]:
    order.append( int(o) )

f = open(sys.argv[1])
f.readline()
for line in f.xreadlines():
    sp = map(removeQuotes, line.split('\t'))
    line = ''
    for o in order:
        line = line + sp[o] + '\t'
    print line.strip()
f.close()
