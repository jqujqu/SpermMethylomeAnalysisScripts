#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export SCRIPT_DIR=${PROJ_DIR}/scripts
export LC_ALL=C

###!!! specify path to TSS file!!!
TSS=~/Annotation/EnsemblGenes_hg19.TSS.bed.unique

mkdir -p ${PROJ_DIR}/hg19Multiz7orthMethylomes/events/filtered/dist_to_TSS
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes/events/filtered/dist_to_TSS

for i in `ls ../*.birth ../*.extension  ../*.contraction ../*.death`; do
name=`basename $i`;
echo $name
bedtools closest -a $i -b ${TSS} -d | \
sort -k1,1 -k2,2g -u -o $name.dist2TSS
done

# plot (Figure 2E)
Rscript ${SCRIPT_DIR}/dist2TSS_violinplot.R




