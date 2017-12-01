#!/bin/bash

# axt alignment files  (https://genome.ucsc.edu/goldenPath/help/axt.html)

ALI=hg19
TAR=mm10
zcat ${ALI}.${TAR}.net.axt.gz | grep chr > ${ALI}.${TAR}.net.axt.bed33 
awk '{OFS="\t"; print $2,$3,$4,"X",0, "+"}' < ${ALI}.${TAR}.net.axt.bed33  | \
LC_ALL=C sort -k1,1 -k3,3g -k2,2g  /dev/stdin  | \
sortbed -c /dev/stdin | \
awk 'BEGIN{OFS="\t";chrom="";start=0; end=0;} 
{if(chrom==$1 && end+1==$2) end=$3; 
 else {if(end-start >0) print chrom,start,end,"X", 0, "+"; 
   chrom=$1; start=$2; end=$3;
 }}
END{ print chrom,start,end,"X", 0, "+" }' < /dev/stdin > ${ALI}.${TAR}.net.bed
