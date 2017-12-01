#!/usr/bin/env python

import sys

genome = dict([(i.strip().split()[0], int(i.strip().split()[-1])) for i in open(sys.argv[1])])
beds = [tuple([i.strip().split()[0], int(i.strip().split()[1]), int(i.strip().split()[2])]) for i in open(sys.argv[2])]

curr = ''
first = True
for i in beds:
    chrom = i[0]
    if chrom != curr:
        if not first:
            print curr, prev, genome[curr]
        first = False
        print chrom, 0, i[1]
    else:
        print curr, prev, i[1]
    curr = chrom
    prev = i[2]
