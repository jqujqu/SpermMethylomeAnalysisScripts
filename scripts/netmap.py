#!/usr/bin/env python
import sys
import os
import re
import argparse

def parse_methcount(line):
  a = line.split()
  if a[2] == "+" :
    return a
  else :
    col4 = a[3].split(":")
    return [a[0], a[1], a[5], col4[0], a[4], col4[1]]


def main():
    
  parser =  argparse.ArgumentParser(description='Map cpg according to netmap')
  parser.add_argument("-n", "--netmap", 
                      help="name of netmap file: e.g. mm9_cpg.hg19net.map")
  parser.add_argument("-o", "--outfile", help="name of output bed file")
  parser.add_argument('meth', help='input methcount file: e.g. mm9.meth')
  args = parser.parse_args()

  out = open(args.outfile, 'w')
  meth = open(args.meth, 'r')
  netmap = open(args.netmap, 'r')

  methline = meth.readline()
  netline = netmap.readline()

  while(len(methline) and len(netline) ):
    [s_chr, s_start, s_strand, s_txt, m, cov] =  parse_methcount(methline)
    [schr, sstart, send, sstrand, tchr, tstart, tend, tstrand] = netline.split() 
    if(s_chr == schr and s_start == sstart):
      out.write(tchr + '\t' + tstart + '\t' + tstrand + '\t' + s_txt +  \
                '\t' + m + '\t' + cov + '\n')
      methline = meth.readline()
      netline = netmap.readline()
    else:
      while(len(methline) and s_chr == schr and int(s_start) < int(sstart)):
        methline = meth.readline()
        if(len(methline)):
          [s_chr, s_start, s_strand, s_txt, m, cov] = parse_methcount(methline)
      while(len(netline) and s_chr == schr and int(s_start) > int(sstart)): 
        netline = netmap.readline()
        if(len(netline)):
          [schr, sstart, send, sstrand, tchr, tstart, tend, tstrand] = netline.split()
      while(len(methline) and s_chr < schr):
        methline = meth.readline()
        if(len(methline)):
          [s_chr, s_start, s_strand, s_txt, m, cov] = parse_methcount(methline)
      while(len(netline) and s_chr > schr):
        netline = netmap.readline()
        if(len(netline)):
          [schr, sstart, send, sstrand, tchr, tstart, tend, tstrand] = netline.split()

  out.close()
  
if __name__ == "__main__":
  main()   
