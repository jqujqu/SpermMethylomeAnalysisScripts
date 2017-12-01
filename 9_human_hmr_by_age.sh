#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export SCRIPT_DIR=${PROJ_DIR}/scripts
export ADSSRC_DIR=~/adssrc
export RMAP_DIR=~/rmapbs/bin
export LC_ALL=C

mkdir -p ${PROJ_DIR}/hg19Multiz7orthMethylomes/events/Human_HMR_by_age_group/
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes/events/Human_HMR_by_age_group/

Human_Meth=${PROJ_DIR}/hg19Multiz7orthMethylomes/Human_Sperm_hg19.meth
Human_HMR=../7sp.Human.HYPO
orth=${PROJ_DIR}/hg19Multiz7orthMethylomes/7orth.bed
desert=${PROJ_DIR}/hg19Multiz7orthMethylomes/7orth_complement.bed

###### entire HMRs #########
${RMAP_DIR}/mapsifter ${Human_HMR} -t ../7sp.Boreoeutheria.HYPO | \
awk '$3-$2>100{print}' > Human_HYPO_since_Boreoeutheria

${RMAP_DIR}/mapsifter ${Human_HMR} -e -t ../7sp.Boreoeutheria.HYPO  | \
${RMAP_DIR}/mapsifter -t ../7sp.Euarchontoglires.HYPO  /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Dog.HYPO /dev/stdin | \
awk '$3-$2>100{print}' > Human_HYPO_since_Euarchontoglires

${RMAP_DIR}/mapsifter ${Human_HMR} -e -t ../7sp.Boreoeutheria.HYPO  | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Euarchontoglires.HYPO  /dev/stdin | \
${RMAP_DIR}/mapsifter -t ../7sp.Catarrhini.HYPO  /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Murinae.HYPO /dev/stdin | \
awk '$3-$2>100{print}' > Human_HYPO_since_Catarrhini

${RMAP_DIR}/mapsifter ${Human_HMR} -e -t ../7sp.Boreoeutheria.HYPO  | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Euarchontoglires.HYPO  /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Catarrhini.HYPO  /dev/stdin | \
${RMAP_DIR}/mapsifter -t ../7sp.Homininae.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -t ../7sp.Gorilla.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Rhesus.HYPO /dev/stdin | \
awk '$3-$2>100{print}' > Human_HYPO_since_Homininae

${RMAP_DIR}/mapsifter ${Human_HMR} -e -t ../7sp.Boreoeutheria.HYPO | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Euarchontoglires.HYPO  /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Catarrhini.HYPO  /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Homininae.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -t ../7sp.Hominini.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -t ../7sp.Chimp.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Gorilla.HYPO /dev/stdin | \
awk '$3-$2>100{print}' > Human_HYPO_since_Hominini

${RMAP_DIR}/mapsifter ${Human_HMR} -e -t ../7sp.Boreoeutheria.HYPO | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Catarrhini.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Homininae.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Hominini.HYPO  /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Chimp.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Gorilla.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Rhesus.HYPO /dev/stdin | \
${RMAP_DIR}/mapsifter -e -t ../7sp.Murinae.HYPO /dev/stdin | \
awk '$3-$2>100{print}' > Human_HYPO_since_Human

for type in `echo Boreoeutheria Euarchontoglires Catarrhini Homininae Hominini Human`; do
echo $type
${SCRIPT_DIR}/fa_interval_CpG_obex.py -b Human_HYPO_since_${type} \
-d /home/rcf-40/jqu/cmb-01/cache/human/hg19/chroms/ \
-o Human_HYPO_since_${type}.cpgobex
sort -k1,1 -k2,2g Human_HYPO_since_${type}.cpgobex -o obextmp && \
mv obextmp  Human_HYPO_since_${type}.cpgobex
roimethstat -L -M Human_HYPO_since_${type} ${Human_Meth} \
-o Human_HYPO_since_${type}.human_sperm_roi
${RMAP_DIR}/mapsifter Human_HYPO_since_${type}.cpgobex \
-t Human_HYPO_since_${type}.human_sperm_roi  \
-o Human_HYPO_since_${type}.cpgobex.filtered && \
mv Human_HYPO_since_${type}.cpgobex.filtered Human_HYPO_since_${type}.cpgobex
${RMAP_DIR}/mapsifter Human_HYPO_since_${type}.human_sperm_roi \
-t Human_HYPO_since_${type}.cpgobex.filtered \
-o Human_HYPO_since_${type}.human_sperm_roi.filtered && \
mv Human_HYPO_since_${type}.human_sperm_roi.filtered Human_HYPO_since_${type}.human_sperm_roi
done

# plot (Figure 4A, S7A)
Rscript ${SCRIPT_DIR}/plot_cpgdensity_HMRage.R

