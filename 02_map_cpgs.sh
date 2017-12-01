#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export METHPIPE_DIR=~/methpipe/bin
export SCRIPT_DIR=${PROJ_DIR}/scripts

## methylome aligned to hg19
#WD=${PROJ_DIR}/hg19MultizMethylomes
#CPGMAP_DIR=${PROJ_DIR}/hg19Multiz100wayMaf/cpgmaps

## methylome aligned to hg19 in 7-way orthologous region
WD=${PROJ_DIR}/hg19Multiz7orthMethylomes
CPGMAP_DIR=${PROJ_DIR}/hg19Multiz100wayMaf/cpgmaps_7orth

SPECIES_LIST=(Human Chimp Gorilla Rhesus Mouse Rat Dog)
FROM_LIST=(hg19 panTro4 gorGor3 rheMac3 mm10 rn5 canFam3)

TO=hg19
export LC_ALL=C

cd $WD
for i in `seq 0 6`; do
SPECIES=${SPECIES_LIST[$i]}
FROM=${FROM_LIST[$i]}
MAP=${CPGMAP_DIR}/${FROM}_cpg_map_to_${TO}.cpgmap
INPUT_DIR=${PROJ_DIR}/OriginalMethylomes/${SPECIES}_${FROM}
SAMPLE=${SPECIES}_Sperm

# collect methylation state, coverage, level and posterior probabilities
${METHPIPE_DIR}/hmr -o ${SAMPLE}.hmr -post-meth ${SAMPLE}.methposterior ${SAMPLE}.meth
bedtools intersect -a ${SAMPLE}.methposterior -b ${SAMPLE}.hmr -wao | \
awk 'BEGIN{OFS="\t"}
     {R=$10; if (R==".") R="HYPER"; split($4,a,":");
      print $1,$2,$3,R"_"a[1]"_"a[2]"_"a[3], $5,$6}
    '> ${SAMPLE}.state_level_post

# map CpGs to new assembly
${SCRIPT_DIR}/netmap.py -n ${MAP} -o ${SAMPLE}_${TO}.tmp ${INPUT_DIR}/${SAMPLE}.state_level_post
grep -v "CpG_0_0" ${SAMPLE}_${TO}.tmp | \
sort -k1,1 -k2,2g -k3,3 /dev/stdin -o ${SAMPLE}_${TO}.state_level_post \
&& rm ${SAMPLE}_${TO}.tmp

# extract HMRs
bash ${SCRIPT_DIR}/extractHMR.sh ${SAMPLE}_${TO}.state_level_post ${SAMPLE}_${TO}.hmr

# extract methcount
awk 'BEGIN{OFS="\t"}
{split($4,a,"_"); print $1,$2,$3, "CpG", $5,a[3]+a[4]}
' < ${SAMPLE}_${TO}.state_level_post > ${SAMPLE}_${TO}.meth

# extract methposterior into bed file
awk '
BEGIN{OFS="\t"}
{split($4, a, "_")}
{print $1,$2,$2+1, a[2]":"a[3]":"a[4], $6, $3}
' < ${SAMPLE}_${TO}.state_level_post > ${SAMPLE}_${TO}.methposterior
done

