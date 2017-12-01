#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export ADSSRC_DIR=~/adssrc
export RMAP_DIR=~/rmapbs/bin
export SCRIPT_DIR=${PROJ_DIR}/scripts
export LC_ALL=C
mkdir -p ${PROJ_DIR}/LineageMutation/regions
mkdir -p ${PROJ_DIR}/LineageMutation/regions_mutation

MutDIR=${PROJ_DIR}/LineageMutation/processed_data
REGDIR=${PROJ_DIR}/LineageMutation/regions

##################################################
# Define regions
cd ${PROJ_DIR}/LineageMutation/regions

for i in `ls ${PROJ_DIR}/hg19Multiz7orthMethylomes/*_Sperm_hg19.hmr`; do
bedtools intersect -a $i \
-b ${PROJ_DIR}/hg19Multiz7orthMethylomes/7orth.bed | \
sort -k1,1 -k2,2g -k3,3g > $(basename $i)
done

${ADSSRC_DIR}/collapsebed *Sperm_hg19.hmr -f -o Sperm.hmr.collapsebed
awk '$5==1{print }' Sperm.hmr.collapsebed | sort -k1,1 -k2,2g -k3,3g > tmp_1x
for sp in `echo Human Chimp Gorilla Rhesus Mouse Rat`; do
bedtools closest -a tmp_1x -b ${sp}_Sperm_hg19.hmr -d | \
awk '$13==0 && ($3-$2)/($9-$8)>0.8{print }' | cut -f 7-12 | \
sort -k1,1 -k2,2g -k3,3g -u > ${sp}_Sperm.hmr.specific
done && rm tmp_1x

# lineage specific HMR
${ADSSRC_DIR}/collapsebed -f -c 4 -o Primate_core.hmr  \
Human_Sperm_hg19.hmr Chimp_Sperm_hg19.hmr Gorilla_Sperm_hg19.hmr Rhesus_Sperm_hg19.hmr
${ADSSRC_DIR}/collapsebed -f -c 2 -o Rodent_core.hmr \
Mouse_Sperm_hg19.hmr Rat_Sperm_hg19.hmr

${RMAP_DIR}/mapsifter  Primate_core.hmr -e -t  Mouse_Sperm_hg19.hmr  | \
${RMAP_DIR}/mapsifter -e -t Rat_Sperm_hg19.hmr /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t Dog_Sperm_hg19.hmr /dev/stdin -o Primate_core.hmr.specific

${RMAP_DIR}/mapsifter  Rodent_core.hmr -e -t  Human_Sperm_hg19.hmr  | \
${RMAP_DIR}/mapsifter -e -t Chimp_Sperm_hg19.hmr /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t Gorilla_Sperm_hg19.hmr /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t Rhesus_Sperm_hg19.hmr /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t Dog_Sperm_hg19.hmr /dev/stdin -o Rodent_core.hmr.specific

# species specific non HMR
for sp in `echo Human Chimp Gorilla Rhesus Mouse Rat`;
do
awk '$5==6{print }' Sperm.hmr.collapsebed | \
${RMAP_DIR}/mapsifter -e -t ${sp}_Sperm_hg19.hmr /dev/stdin -o ${sp}_Sperm.nonhmr.specific;
done

${ADSSRC_DIR}/collapsebed -c 3 -o MRD.bed  \
Mouse_Sperm_hg19.hmr Rat_Sperm_hg19.hmr Dog_Sperm_hg19.hmr
${ADSSRC_DIR}/collapsebed -o Primate.collapsed.bed \
Human_Sperm_hg19.hmr Chimp_Sperm_hg19.hmr Gorilla_Sperm_hg19.hmr Rhesus_Sperm_hg19.hmr
${RMAP_DIR}/mapsifter MRD.bed -e -t Primate.collapsed.bed -o Primate_core.nonhmr.specific

${ADSSRC_DIR}/collapsebed -c 5 -o HCGRD.bed  \
Human_Sperm_hg19.hmr Chimp_Sperm_hg19.hmr Gorilla_Sperm_hg19.hmr \
Rhesus_Sperm_hg19.hmr Dog_Sperm_hg19.hmr

${ADSSRC_DIR}/collapsebed -o Rodent.collapsed.bed Mouse_Sperm_hg19.hmr Rat_Sperm_hg19.hmr
${RMAP_DIR}/mapsifter HCGRD.bed -e -t Rodent.collapsed.bed -o Rodent_core.nonhmr.specific

rm MRD.bed Primate.collapsed.bed HCGRD.bed Rodent.collapsed.bed

##################################################
# count mutations
cd ${PROJ_DIR}/LineageMutation/regions_mutation

cat < node_num_species.txt <<-END
0       Human
1       Hominini
2       Chimp
3       Homininae
4       Gorilla
5       Catarrhini
6       Rhesus
7       Euarchontoglires
8       Mouse
9       Murinae
10      Rat
11      Boreoeutheria
12      Dog
END

# region files for lineage-specifc HMRs
echo "Region Total_size Aligned_size"| \
tr ' ' '\011' < /dev/stdin > region_size.txt

for i in `ls ${REGDIR}/*.hmr.specific`; do
tot_size=`awk '{i+=$3-$2} END{print i/1e6}' < $i`
aligned_size=`bedtools intersect -sorted -a $i \
-b ${MutDIR}/Compara.15_eutherian_mammals_EPO.aligned.bed | \
awk '{i+=$3-$2} END{print i/1e6}' `;
echo $(basename $i) $tot_size $aligned_size | \
tr ' ' '\011' < /dev/stdin
done >> region_size.txt

for i in `ls ${REGDIR}/*.hmr.specific `; do
bash ${SCRIPT_DIR}/CallMut.sh $i  $(basename $i).txt
done

for i in `ls ${REGDIR}/*.nonhmr.specific `; do
bash ${SCRIPT_DIR}/CallMut.sh $i  $(basename $i).txt
done

# plot (Figure 4B)
Rscript ${SCRIPT_DIR}/LineageMutationPlot.R
