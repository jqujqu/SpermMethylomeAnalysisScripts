#!/bin/bash
awk 'BEGIN{prev=""; d = 1000; OFS="\t"}
 {split($4,a,"_");
  if (substr(a[1],1,4)=="HYPO" && prev!= a[1]) {
    if (substr(prev,1,4)=="HYPO" && count >=5) {
      print chrom, start, end, prev, count, "."
    }
    chrom=$1; start=$2; end=$2+1; prev=a[1]; count=1;
  } else if (substr(a[1],1,4)=="HYPO" && prev== a[1] && $2-end <= d && chrom==$1) {
    count += 1; end=$2+1;
  } else if (substr(a[1],1,4)=="HYPO" && prev== a[1] && ($2-end > d || chrom!=$1)) {
    if (count>=5) print chrom, start, end, prev, count, ".";
    chrom=$1; start=$2; end=$2+1; prev=a[1]; count=1;
  } else if (substr(a[1],1,4)!="HYPO" && substr(prev,1,4)== "HYPO") {
    if (count>=5) print chrom, start, end, prev, count, ".";
    prev = a[1]
  }
}
END{ if (substr(prev,1,4)== "HYPO" && count>=5) print chrom, start, end, prev, count, ".";}
' < $1  > $2
