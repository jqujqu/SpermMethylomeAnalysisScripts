#!/usr/bin/env python
import sys
import os
import re

def read_bed_interval(bedfile):
  f = open(bedfile, 'r')
  intervals = {}
  for line in f :
    [chrom, start, end ] = line.split()[:3]
    if chrom in intervals :
      intervals[chrom].append([int(start), int(end)])
    else:
      intervals[chrom]=[[int(start), int(end)]]
  return intervals


####   obex  ########
def obex(name, seq, interval, out):
  start = interval[0]
  end = interval[1]
  s = seq[start:end].upper()
  gc = float(s.count('C')+s.count('G'))
  CpG = s.count('CG')
  if gc > 0:
    cpg = float(CpG)
    effective_size = len(s)-s.count('N')
    gc = gc/(effective_size)/2
    expcpg = gc*gc*len(s)
    obexval = cpg/expcpg
    CorG = s.count('C')+s.count('G')
    out.write(name + "\t" + str(start)+"\t" + str(end) + "\tCpG:" + 
              str(CpG) + ":"+ str(CorG) + ":" + str(effective_size) + 
               "\t" + str(obexval) + "\t+\n")
  return

def obex_mask(name, seq, interval, out):
  start = interval[0]
  end = interval[1]
  s = seq[start:end]
  gc = float(s.count('C')+s.count('G'))
  CpG = s.count('CG')
  if gc > 0:
    cpg = float(CpG)
    effective_size = float(s.count('A')+s.count('T')) + gc
    gc = gc/(effective_size)/2
    expcpg = gc*gc*len(s)
    obexval = cpg/expcpg
    CorG = s.count('C')+s.count('G')
    out.write(name + "\t" + str(start)+"\t" + str(end) + "\tCpG:" +
              str(CpG) + ":"+ str(CorG) + ":" + str(effective_size) +
               "\t" + str(obexval) + "\t+\n")
  return


def main(argv):
  import argparse
  parser = argparse.ArgumentParser(description='Compute CpG ob/ex ratio in genomic intervals')
  parser.add_argument('-b', metavar='BED', dest='interval',action='store', help='Interval BED file')
  parser.add_argument('-d', default="./", metavar='PATH', dest='fapath', help='Directory containing fa files')
  parser.add_argument('--mask', dest='mask', action='store_const',
                   const=True, default=False,
                   help='only use unmasked sequence (default: use all sequence)')
  parser.add_argument('-o', metavar='BED', dest = 'obexfile', action='store',
                      help='Output bed file Name=#CpG:#CorG')
  args = parser.parse_args()
                      
  intervalfile = args.interval
  mask = args.mask
  intervals = read_bed_interval(intervalfile)
  fapath = args.fapath
  outfile = args.obexfile
  out = open(outfile, 'w') 
  
  from os import listdir
  from os.path import isfile, join
  fafiles = [ f for f in listdir(fapath) if isfile(join(fapath,f)) and f.split('.')[-1] == 'fa' ]

  for chrom in intervals :
    fafile=chrom+'.fa'
    assert(fafile in fafiles)
    f = open(join(fapath,fafile), 'rU')
    name = chrom
    for line in f :
      if (line[0] == '>'):
        seq = ''
      else :
        seq += line.rstrip('\n')
    for interval in intervals[chrom] :
      if mask:
        obex_mask(name, seq, interval , out)
      else :
        obex(name, seq, interval , out)
  return


if __name__ == "__main__":
  main(sys.argv)

