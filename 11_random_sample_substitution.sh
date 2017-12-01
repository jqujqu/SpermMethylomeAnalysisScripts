#!/bin/bash
export PROJ_DIR=~/sperm_methylome_evolution
export SCRIPT_DIR=${PROJ_DIR}/scripts
export LC_ALL=C
mkdir -p ${PROJ_DIR}/LineageMutation/random_sample_region_0001
cd ${PROJ_DIR}/LineageMutation/random_sample_region_0001
MutDatDIR=${PROJ_DIR}/LineageMutation/Ensembl75_Compara15/processed_data

# randomly sample 0.01% of aligned sequences 5000 times
LINE=`awk 'END{print int(NR/10000)}' < ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.aligned.bed`

for j in `seq 1 5000`; do
shuf -n ${LINE} ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.aligned.bed | \
sort -k1,1 -k2,2g | awk 'BEGIN{OFS="\t";} $3>$2 {print $1,$2,$3,"X",0,"+"}' > tmp
Primate=`bedtools intersect -sorted -a tmp -b ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.mutation.bed.node_5 | wc -l`
Rodent=`bedtools intersect -sorted -a tmp -b ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.mutation.bed.node_9 | wc -l`
Mouse=`bedtools intersect -sorted -a tmp -b ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.mutation.bed.node_8 | wc -l`
Rat=`bedtools intersect -sorted -a tmp -b ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.mutation.bed.node_10 | wc -l`
Human=`bedtools intersect -sorted -a tmp -b ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.mutation.bed.node_0 | wc -l`
Chimp=`bedtools intersect -sorted -a tmp -b ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.mutation.bed.node_2 | wc -l`
Rhesus=`bedtools intersect -sorted -a tmp -b ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.mutation.bed.node_6 | wc -l`
Gorilla=`bedtools intersect -sorted -a tmp -b ${MutDatDIR}/Compara.15_eutherian_mammals_EPO.mutation.bed.node_4 | wc -l`
echo $j $Primate $Rodent $Mouse $Rat $Human $Chimp $Rhesus $Gorilla >> sampling_results.txt
done
