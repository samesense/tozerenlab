#!/usr/bin/python

import sys

if len(sys.argv) == 1:
    print 'enter a list of lists'
    sys.exit(0)

intersection = dict()
lsts = []

for fname in sys.argv[1:]:
    f = open(fname)
    lsts.append(dict())
    for line in f.xreadlines():
        lsts[-1][line.strip()] = True
    f.close()

for k in lsts[0].keys():
    add = True
    for d2 in lsts[1:]:
        if not d2.has_key(k):
            add = False
            break
    if add:
        print k
