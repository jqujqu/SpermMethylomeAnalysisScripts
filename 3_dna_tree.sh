#!/bin/bash
export PROJ_DIR=/home/rcf-40/jqu/staging/test
export SCRIPT_DIR=${PROJ_DIR}/scripts
export GALAXY_DIR=~/cmb-01/tools/galaxy

mkdir -p ${PROJ_DIR}/DNATree
cd ${PROJ_DIR}/DNATree
for i in `seq 1 22`; do
${GALAXY_DIR}/tools/maf/maf_to_fasta_concat.py \
hg19,panTro4,gorGor3,rheMac3,mm10,rn5,canFam3 \
${PROJ_DIR}/hg19Multiz100wayMaf/cpgmaps/chr${i}.7sp.maf chr${i}.7sp.fasta
done

# estimated unrestricted tree by chromosome
Rscript ${SCRIPT_DIR}/rphast.R

