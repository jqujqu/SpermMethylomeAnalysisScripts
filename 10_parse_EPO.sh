#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export SCRIPT_DIR=${PROJ_DIR}/scripts
export LC_ALL=C

mkdir -p ${PROJ_DIR}/LineageMutation/processed_data
cd ${PROJ_DIR}/LineageMutation/processed_data

for i in `ls ../raw_data/Compara.15_eutherian_mammals_EPO.chr*_*.emf.gz`; do
f=$(basename $i)
# target species hard-coded
# node numbered in
${SCRIPT_DIR}/parse_emf_readall.py -o ${f/.emf.gz/} $i
sort -k1,1 -k2,2g -k3,3g ${f/.emf.gz/}.mutation.bed -o ${f/.emf.gz/}.sorted && \
mv ${f/.emf.gz/}.sorted ${f/.emf.gz/}.mutation.bed
sort -k1,1 -k2,2g -k3,3g ${f/.emf.gz/}.aligned.bed  -o ${f/.emf.gz/}.sorted && \
mv ${f/.emf.gz/}.sorted ${f/.emf.gz/}.aligned.bed
done

sort -k1,1 -k2,2g -k3,3g -m Compara.15_eutherian_mammals_EPO.chr*_*.mutation.bed  \
-o Compara.15_eutherian_mammals_EPO.mutation.bed
sort -k1,1 -k2,2g -k3,3g -m Compara.15_eutherian_mammals_EPO.chr*_*.aligned.bed \
-o Compara.15_eutherian_mammals_EPO.aligned.bed

awk '{split($4,a,":");
      print > "Compara.15_eutherian_mammals_EPO.mutation.bed."a[1]}
' < Compara.15_eutherian_mammals_EPO.mutation.bed



