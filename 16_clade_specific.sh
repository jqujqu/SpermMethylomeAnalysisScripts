#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export SCRIPT_DIR=${PROJ_DIR}/scripts
export ADSSRC_DIR=~/adssrc
export RMAP_DIR=~/rmapbs/bin
export LC_ALL=C

ORTH_DIR=${PROJ_DIR}/hg19Multiz7orthMethylomes
CLADE_DIR=${PROJ_DIR}/CladeSpecific

mkdir -p $ORTH_DIR $CLADE_DIR

# specify path to annotation file
Promoter=~/Annotation/EnsemblGenes_hg19.TSSpm1k_merged.bed

#######
cd ${ORTH_DIR}
${ADSSRC_DIR}/collapsebed *Sperm_hg19.hmr -f -c 7 -o hg19_sperm_HMR_shared_Boreoeutheria.bed

${ADSSRC_DIR}/collapsebed -f -c 6 *Sperm_hg19.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Dog_Sperm_hg19.hmr -o hg19_sperm_HMR_shared_Euarchontoglire.bed

${ADSSRC_DIR}/collapsebed -f -c 4 *Sperm_hg19.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Dog_Sperm_hg19.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Mouse_Sperm_hg19.hmr |\
${RMAP_DIR}/mapsifter /dev/stdin -e -t Rat_Sperm_hg19.hmr  \
-o hg19_sperm_HMR_shared_Catarrhini.bed

${ADSSRC_DIR}/collapsebed -f -c 3 *Sperm_hg19.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -t Human_Sperm_hg19.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -t Chimp_Sperm_hg19.hmr |\
${RMAP_DIR}/mapsifter /dev/stdin -t Gorilla_Sperm_hg19.hmr  \
-o hg19_sperm_HMR_shared_Homininae.bed

${ADSSRC_DIR}/collapsebed -f -c 2 *Sperm_hg19.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -t Human_Sperm_hg19.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -t Chimp_Sperm_hg19.hmr \
-o hg19_sperm_HMR_shared_Hominini.bed

${ADSSRC_DIR}/collapsebed -f -c 2 *Sperm*.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -t Mouse_Sperm_hg19.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -t Rat_Sperm_hg19.hmr \
-o hg19_sperm_HMR_shared_Murinae.bed

#######
cd  ${CLADE_DIR}
#${ORTH_DIR}/hg19_sperm_HMR_shared_Boreoeutheria.bed
#${ORTH_DIR}/hg19_sperm_HMR_shared_Euarchontoglire.bed
#${ORTH_DIR}/hg19_sperm_HMR_shared_Catarrhini.bed
#${ORTH_DIR}/hg19_sperm_HMR_shared_Homininae.bed
#${ORTH_DIR}/hg19_sperm_HMR_shared_Hominini.bed
#${ORTH_DIR}/hg19_sperm_HMR_shared_Murinae.bed

hHMR=${ORTH_DIR}/Human_Sperm_hg19.hmr
${RMAP_DIR}/mapsifter $hHMR -t ${ORTH_DIR}/hg19_sperm_HMR_shared_Boreoeutheria.bed -o Human_Sperm_HMR.since_Boreoeutheria
${RMAP_DIR}/mapsifter $hHMR -t ${ORTH_DIR}/hg19_sperm_HMR_shared_Euarchontoglire.bed | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Boreoeutheria -o Human_Sperm_HMR.since_Euarchontoglire
${RMAP_DIR}/mapsifter $hHMR -t ${ORTH_DIR}/hg19_sperm_HMR_shared_Catarrhini.bed | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Boreoeutheria  | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Euarchontoglire -o Human_Sperm_HMR.since_Catarrhini
${RMAP_DIR}/mapsifter $hHMR -t ${ORTH_DIR}/hg19_sperm_HMR_shared_Homininae.bed |\
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Catarrhini | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Boreoeutheria  | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Euarchontoglire -o Human_Sperm_HMR.since_Homininae
${RMAP_DIR}/mapsifter $hHMR -t ${ORTH_DIR}/hg19_sperm_HMR_shared_Hominini.bed | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Homininae |\
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Catarrhini | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Boreoeutheria  | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Euarchontoglire -o Human_Sperm_HMR.since_Hominini
${RMAP_DIR}/mapsifter  $hHMR -t ${ORTH_DIR}/Human_Sperm_hg19.hmr | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Hominini | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Homininae |\
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Catarrhini | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Boreoeutheria  | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Human_Sperm_HMR.since_Euarchontoglire -o Human_Sperm_HMR.since_Human

mHMR=${ORTH_DIR}/Mouse_Sperm_hg19.hmr
${RMAP_DIR}/mapsifter $mHMR -t ${ORTH_DIR}/hg19_sperm_HMR_shared_Boreoeutheria.bed -o Mouse_Sperm_HMR.since_Boreoeutheria
${RMAP_DIR}/mapsifter $mHMR -t ${ORTH_DIR}/hg19_sperm_HMR_shared_Euarchontoglire.bed | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Mouse_Sperm_HMR.since_Boreoeutheria -o Mouse_Sperm_HMR.since_Euarchontoglire
${RMAP_DIR}/mapsifter $mHMR -t ${ORTH_DIR}/hg19_sperm_HMR_shared_Murinae.bed | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Mouse_Sperm_HMR.since_Euarchontoglire | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Mouse_Sperm_HMR.since_Boreoeutheria -o Mouse_Sperm_HMR.since_Murinae
${RMAP_DIR}/mapsifter $mHMR -e -t ${ORTH_DIR}/hg19_sperm_HMR_shared_Murinae.bed | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Mouse_Sperm_HMR.since_Euarchontoglire | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t Mouse_Sperm_HMR.since_Boreoeutheria -o Mouse_Sperm_HMR.since_Mouse

desert=${ORTH_DIR}/7orth_complement.bed
for i in `echo Boreoeutheria Euarchontoglire Catarrhini Homininae Hominini Human`; do
${RMAP_DIR}/mapsifter Human_Sperm_HMR.since_${i}  -e -t $desert -o tmp && mv tmp Human_Sperm_HMR.since_${i}
done
for i in `echo Boreoeutheria Euarchontoglire Murinae Mouse`; do
${RMAP_DIR}/mapsifter Mouse_Sperm_HMR.since_${i} -e -t $desert -o tmp && mv tmp Mouse_Sperm_HMR.since_${i}
done


for i in `echo Boreoeutheria Euarchontoglire Catarrhini Homininae Hominini Human`; do
${RMAP_DIR}/mapsifter Human_Sperm_HMR.since_${i} -t $Promoter \
-o Human_Sperm_HMR.since_${i}.promoter
${RMAP_DIR}/mapsifter Human_Sperm_HMR.since_${i} -e -t $Promoter \
-o Human_Sperm_HMR.since_${i}.nonpromoter
done
for i in `echo Boreoeutheria Euarchontoglire Murinae Mouse`; do
${RMAP_DIR}/mapsifter Mouse_Sperm_HMR.since_${i} -t $Promoter  \
-o Mouse_Sperm_HMR.since_${i}.promoter
${RMAP_DIR}/mapsifter Mouse_Sperm_HMR.since_${i} -e -t $Promoter \
-o Mouse_Sperm_HMR.since_${i}.nonpromoter
done

##################################################
# plot (Figure 3A)
Rscript ${PROJ_DIR}/scripts/clade_specific_plot.R
