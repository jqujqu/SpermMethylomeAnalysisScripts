#!/usr/bin/python

import gzip
import sys
import argparse
import numpy as np
import re

###############################################################################
"""
maf parsing
"""

def isGap(s):
  return s[0] == '-'

def isBlockStart(s):
  return s.split()[0] == 'a'

def isSequence(s) :
  return s.split()[0]=='s'

def readMafBlock(f) :
  line = f.readline()
  while len(line) and not isBlockStart(line):
    line = f.readline()
  score = 0
  for field in line.split() :
    if field[:6] == "score=" :
      score = float(field[6:])
  seqlines =[]
  line = f.readline()
  while (line != '\n') :
    if isSequence(line):
      seqlines.append(line)
    line = f.readline()
  return (seqlines, score)

def countGap(seq) :
  gaps = [1 if x == '-' else 0 for x in seq]
  gap_acc =  np.cumsum(gaps).tolist()
  return gap_acc

class SeqRec:
  """
  s hg16.chr7    27707221 13 + 158545518 gcagctgaaaaca
  """
  def __init__(self, seq):
    s = seq.split()
    self.src = s[1]
    self.start = int(s[2])
    self.size = int(s[3])
    self.strand = s[4]
    self.srcSize = int(s[5])
    self.text = s[6].upper()
    self.gapAcc = countGap(self.text)
    self.srcStart = 0
    if self.strand[0] == '+' :
      self.srcStart = int(s[2])
    else :
      self.srcStart = self.srcSize - self.start -self.size +1

def peek_line(f):
  pos = f.tell()
  line = f.readline()
  f.seek(pos)
  return line

def parse_CpG_pos(alignment, score, output) :
  """
  alignment is a list of SeqRec
  """
  nsp = len(alignment)
  for i in range(nsp) :
    for m in re.finditer(r"C-*G", alignment[i].text) :
      start_in_seq = m.start()
      end_in_seq = m.end() - 1
      Cpos = alignment[i].start + start_in_seq - \
             alignment[i].gapAcc[start_in_seq]
      Gpos = alignment[i].start + end_in_seq - \
             alignment[i].gapAcc[end_in_seq]
      Crefpos = alignment[0].start + start_in_seq - \
                alignment[0].gapAcc[start_in_seq]
      Grefpos = alignment[0].start + end_in_seq - \
                alignment[0].gapAcc[end_in_seq]
      if alignment[i].strand == "-" :
          Cpos = alignment[i].srcSize - Cpos - 1
          Gpos = alignment[i].srcSize - Gpos - 1
      if alignment[0].strand == "-" :
          Crefpos = alignment[0].srcSize - Crefpos -1
          Grefpos = alignment[0].srcSize - Grefpos -1
      if not (isGap(alignment[0].text[start_in_seq]) and
              isGap(alignment[0].text[end_in_seq])) :
        output.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(alignment[i].src,
        min(Cpos,Gpos), max(Cpos, Gpos), m.group(0), alignment[i].strand,
        alignment[0].src, min(Crefpos,Grefpos), max(Crefpos, Grefpos),
        alignment[0].text[m.start():m.end()], alignment[0].strand, score) + '\n')


def pass_filter(alignment, splst) :
  if len(splst) == 0 :
      return True
  nsp = len(alignment)
  sp = set([str.split(seqrec.src, ".")[0] for seqrec in alignment])
  return sp == set(splst)


###############################################################################

def main():
  parser = argparse.ArgumentParser(description='Parse locations of CpG sites in'
                                    'aligning sequence and reference sequence',
                                   prog='mafFindCpG.py')
  parser.add_argument('--maf', dest='maf', help='maf file')
  parser.add_argument("--splst", dest='splst', default='', 
                      help='species list, require block to have one ' 
                      'sequence record of each species (optional)')
  parser.add_argument('--output', dest='output', help='output file')
  if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)
  args = parser.parse_args()
  maf = args.maf  #"chr1.7sp.maf"
  out = open(args.output, 'w')  # "XXX.CpG.tsv"

  splst = []
  if len(args.splst) :
    splst = [line.rstrip('\n') for line in open(args.splst)]

  f = gzip.open(maf, 'r') if maf[-2:]=="gz" else open(maf, 'r')
  line =  peek_line(f)
  alignments = []
  while len(line):
    alignment =[]
    (seq_lines, score) = readMafBlock(f)
    for seq in seq_lines :
      alignment.append(SeqRec(seq))
      if pass_filter(alignment, splst) :
        parse_CpG_pos(alignment, score, out)
    line = peek_line(f)
  f.close()
  out.close()


if __name__ == "__main__":
  main()

