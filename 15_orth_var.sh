#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export RMAP_DIR=~/rmapbs/bin
export SCRIPT_DIR=${PROJ_DIR}/scripts
export LC_ALL=C

####### for each species (e.g. Chimp panTro4)#############
cd  ${PROJ_DIR}/OriginalOrthVar/OrthVarChimp_panTro4
HMR=${PROJ_DIR}/OriginalMethylomes/Chimp_panTro4/Chimp_Sperm.hmr

ASS=panTro4
cpgmap=${PROJ_DIR}/hg19Multiz100wayMaf/cpgmaps_7orth/${ASS}_cpg_map_to_hg19.cpgmap
HMR=${PROJ_DIR}/OriginalMethylomes/Chimp_panTro4/Chimp_Sperm.hmr
bedtools intersect -u -a $HMR  -b $cpgmap > $(basename $HMR).ORTH
${RMAP_DIR}/mapsifter -e -t $(basename $HMR).ORTH $HMR -o $(basename $HMR).VAR

# specify paths to annotation files
Repeat=~/panfs/cache/chimp/panTro4/Annotation/RepeatMasker/RepeatMasker_panTro4.bed.merged
Promoter=~/panfs/cache/chimp/panTro4/Annotation/EnsemblGenes/EnsemblGenes_panTro4.TSS.bed.pm1k.merged
for i in ` echo ORTH VAR`; do
${RMAP_DIR}/mapsifter $(basename $HMR).${i} -t $Promoter -o $(basename $HMR).${i}_promoter
${RMAP_DIR}/mapsifter $(basename $HMR).${i} -e -t $Promoter | \
${RMAP_DIR}/mapsifter /dev/stdin -t ${Repeat} -o $(basename $HMR).${i}_repeat
${RMAP_DIR}/mapsifter $(basename $HMR).${i} -e -t $Promoter | \
${RMAP_DIR}/mapsifter /dev/stdin -e -t ${Repeat} -o $(basename $HMR).${i}_other
done

##################################################
# plot (Figure 3B)
Rscript ${PROJ_DIR}/scripts/orth_var_plot.R
