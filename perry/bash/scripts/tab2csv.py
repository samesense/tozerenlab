#!/usr/bin/python
import sys

for line in sys.stdin.xreadlines():
    print line.replace('\t',',').strip()
