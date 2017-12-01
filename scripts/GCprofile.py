#!/usr/bin/env python

import sys
import os
import re
import argparse

def read_bed_interval(bedfile):
  f = open(bedfile, 'r')
  intervals = {}
  for line in f :
    [chrom, start, end, name, center, strand ] = line.split()
    strand = 2*int(strand=="+")-1 
    if chrom in intervals :
      intervals[chrom].append([int(start), int(end), int(center), strand])
    else:
      intervals[chrom]=[[int(start), int(end),  int(center), strand]]
  return intervals


def slidingWindow(sequence,winSize,step=1):
  """Returns a generator that will iterate through
  the defined chunks of input sequence.  Input sequence
  must be iterable."""
  # Verify the inputs
  try: it = iter(sequence)
  except TypeError:
    raise Exception("**ERROR** sequence must be iterable.")
  if not ((type(winSize) == type(0)) and (type(step) == type(0))):
    raise Exception("**ERROR** type(winSize) and type(step) must be int.")
  if step > winSize:
    raise Exception("**ERROR** step must not be larger than winSize.")
  if winSize > len(sequence):
    raise Exception("**ERROR** winSize must not be larger than sequence length.")
  # Pre-compute number of chunks to emit
  numOfChunks = ((len(sequence)-winSize)/step)+1
  # Do the work
  for i in range(0,numOfChunks*step,step):
    yield sequence[i:i+winSize]


def obex(seq, interval, halfwinsize, stepsize, profile):
  """ profile is a dictionary, with keys as distances to reference point """
  [start, end, ref, strand] = interval
  s = seq[(start-halfwinsize):(end+halfwinsize)].upper()
  chunks = slidingWindow(s, 2*halfwinsize, stepsize)
  dist = start - ref
  for window in chunks :
    dist += stepsize
    Gsize = window.count('G')
    Csize = window.count('C')
    CpG = window.count('CG')
    d = dist*strand
    tot = len(window)-window.count('N')
    if d in profile.keys(): 
      if strand == 1 :
        profile[d]["G"] += Gsize 
        profile[d]["C"] += Csize
      else :
        profile[d]["G"] += Csize 
        profile[d]["C"] += Gsize
      profile[d]["CpG"] += CpG
      profile[d]["TOT"] += tot
      profile[d]["Pairs"] += tot-1
    else :
      if strand == 1 :
        profile[d] = {"G": Gsize, "C":Csize, "CpG": CpG, "TOT": tot, "Pairs":tot-1}
      else :
        profile[d] = {"C": Gsize, "G":Csize, "CpG": CpG, "TOT": tot, "Pairs":tot-1}
  return


def main(argv): 
  parser = argparse.ArgumentParser(description='Profile CpG ob/ex ratio, GC content, GC skew')
  parser.add_argument('-b', metavar='BED', dest='interval', action='store', 
                      help='Interval BED file: 6-column, reference point position in 5th column')
  parser.add_argument('-w', metavar='half-window-size', dest='halfwinsize',action='store', help='half window size')
  parser.add_argument('-s', metavar='step-size', dest='stepsize',action='store', help='step size')
  parser.add_argument('-f', default="./", metavar='PATH', dest='fapath', help='Directory containing fa files')
  parser.add_argument('--mask', dest='mask', action='store_const', const=True, default=False,
                      help='only use unmasked sequence (default: use all sequence)')
  parser.add_argument('-o', metavar='output', dest = 'outfile', action='store',
                      help='Output format <DIST> <CpGo/e> <GC content> <GC skew>')
  args = parser.parse_args()
  
  intervalfile = args.interval
  halfwinsize = int(args.halfwinsize)
  stepsize = int(args.stepsize)
  mask = args.mask
  intervals = read_bed_interval(args.interval)
  fapath = args.fapath
  outfile = args.outfile
  out = open(outfile, 'w') 
  
  from os import listdir
  from os.path import isfile, join
  fafiles = [ f for f in listdir(fapath) if isfile(join(fapath,f)) and f.split('.')[-1] == 'fa' ]
  
  profile = {}
  for chrom in intervals :
    fafile = chrom+'.fa'
    assert(fafile in fafiles)
    f = open(join(fapath,fafile), 'rU')
    name = chrom
    for line in f :
      if (line[0] == '>'):
        seq = ''
      else :
        seq += line.rstrip('\n')
    for interval in intervals[chrom] :
      obex(seq, interval, halfwinsize, stepsize, profile)
  
  out.write("Dist" + "\t" + "CpGoe" + "\t" + "GCcontent" + "\t" + "GCskew" + "\t" + "Total"+ "\n")
  for d in sorted(profile.keys()):
    stat = profile[d]    
    gorc = float(stat["C"]+stat["G"])/stat["TOT"]/2
    if gorc > 0 :
      CpGoe = str( float(stat["CpG"])/stat["Pairs"]/pow(gorc,2)  )
    else :
      CpGoe = "NA"
    GCcontent = str(float(stat["C"]+stat["G"])/stat["TOT"])
    if stat["C"] + stat["G"] > 0 :
      GCskew = str(float(stat["G"]-stat["C"])/(stat["C"]+stat["G"]))
    else :
      GCskew = "NA"
    out.write(str(d) + "\t" + CpGoe + "\t" + GCcontent + "\t" + GCskew + "\t" + str(stat["Pairs"])+ "\n")
  return

if __name__ == "__main__":
  main(sys.argv)

