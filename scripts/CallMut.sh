#!/bin/bash

# Call this script with 2 arguments:
# ./CallMut.sh regionfile outfile

reg=$1
output=$2

cd ${PROJ_DIR}/LineageMutation/regions_mutation

MutDIR=${PROJ_DIR}/LineageMutation/processed_data

tot=`bedtools intersect -sorted -b ${MutDIR}/Compara.15_eutherian_mammals_EPO.aligned.bed -a ${reg}| awk '{i+=$3-$2} END{print i}'`
echo "branch name tot AC AG AT CG CT GT" > ${output}
while read line ; do 
i=`echo $line | awk '{print $1}'`
branch=Compara.15_eutherian_mammals_EPO.mutation.bed.node_${i}
sub=`bedtools intersect -a ${MutDIR}/${branch}  -b ${reg} | \
tr ':' '\011'  <  /dev/stdin | \
awk 'BEGIN{map["AC"]="AC"; map["AG"]="AG"; map["AT"]="AT"; map["CG"]="CG"; map["CT"]="CT"; map["GT"]="GT";
           map["TG"]="AC"; map["TC"]="AG"; map["TA"]="AT"; map["GC"]="CG"; map["GA"]="CT"; map["CA"]="GT";
           mut["AC"]=0; mut["AG"]=0; mut["AT"]=0; mut["CG"]=0; mut["CT"]=0; mut["GT"]=0;}
{mut[map[substr($5,2,1)substr($6,2,1)]]+=1} END{print mut["AC"], mut["AG"], mut["AT"], mut["CG"], mut["CT"], mut["GT"]}'`
echo $line $tot $sub >> $output
done < node_num_species.txt 

tr ' ' '\011' < $output > ${output}.tmp && mv ${output}.tmp ${output}
