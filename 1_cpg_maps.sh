#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export GALAXY_DIR=~/cmb-01/tools/galaxy
export UCSC_UTILS_DIR=~/cmb-01/tools/GenomeBrowser
export SCRIPT_DIR=${PROJ_DIR}/scripts


# output directory
MAF_DIR=${PROJ_DIR}/hg19Multiz100wayMaf
CPGMAP_7WAY_DIR=${MAF_DIR}/cpgmaps_7orth
CPGMAP_DIR=${MAF_DIR}/cpgmaps

export LC_ALL=C

for i in `echo $(seq 1 22) M X Y`; do
# filter for 7 species (remove lines representing species not from the species.lst)
${UCSC_UTILS_DIR}/mafSpeciesSubset ${MAF_DIR}/chr${i}.maf.gz \
${MAF_DIR}/species.lst ${MAF_DIR}/chr${i}.7sp.maf
done

#############################
# maps for pairs of species #
#############################

for ass in `echo panTro4 gorGor3 rheMac3 mm10 rn5 canFam3`; do
for i in `echo $(seq 1 22) M X Y`; do
# filter and stitch blocks for 2 species
${GALAXY_DIR}/tools/maf/maf_thread_for_species.py ${MAF_DIR}/chr${i}.7sp.maf \
    ${CPGMAP_DIR}/chr${i}.hg19_${ass}.stitched.maf hg19,${ass}
# extract cpgmaps
${SCRIPT_DIR}/mafFindCpGFilter.py --maf ${CPGMAP_DIR}/chr${i}.hg19_${ass}.stitched.maf \
  --output ${CPGMAP_DIR}/chr${i}.hg19_${ass}.stitched.cpg_map_hg19 && \
  rm ${CPGMAP_DIR}/chr${i}.hg19_${ass}.stitched.maf
done
prefix=${CPGMAP_DIR}/${ass}_cpg_map_to_hg19
cat ${CPGMAP_DIR}/chr*.hg19_${ass}.stitched.cpg_map_hg19 |\
  awk '{split($1, a, ".")} a[1]!="hg19"{print}'> ${prefix}
sort -k6,6 -k7,7g -u ${prefix} -o ${prefix}.sorted_hg19
sort -k1,1 -k2,2g -u ${prefix}.sorted_hg19 -o ${prefix}.sorted_${ass}
awk 'BEGIN{OFS="\t"} {split($1,a, "."); split($6,b,".");}
  {print a[2],$2,$3,$5, b[2], $7,$7+1,$10}
  ' < ${prefix}.sorted_${ass} > ${prefix}.cpgmap && \
rm ${prefix} ${prefix}.sorted_hg19 ${prefix}.sorted_${ass} \
  ${CPGMAP_DIR}/chr*.hg19_${ass}.stitched.cpg_map_hg19
done


#########################
# maps for 7-way blocks #
#########################
# generate map in 7-way orthologous blocks
for i in `echo $(seq 1 22) M X Y`; do
# filter for blocks having all of the passed in species and fuse adjacent blocks
${GALAXY_DIR}/tools/maf/maf_thread_for_species.py ${MAF_DIR}/chr${i}.7sp.maf \
  ${CPGMAP_7WAY_DIR}/chr${i}.7orth.stitched.maf hg19,panTro4,gorGor3,rheMac3,mm10,rn5,canFam3
# all-ref
${SCRIPT_DIR}/mafFindCpGFilter.py \
  --splst ${CPGMAP_7WAY_DIR}/species.lst \
  --maf ${CPGMAP_7WAY_DIR}/chr${i}.7orth.stitched.maf \
  --output ${CPGMAP_7WAY_DIR}/chr${i}.7orth.stitched.cpg_map_hg19
done

# merge, sort and filter
cat ${CPGMAP_7WAY_DIR}/chr*.7orth.stitched.cpg_map_hg19 | \
awk -v odir=${CPGMAP_7WAY_DIR} '{split($1,a, ".");
  print > odir"/"a[1]"_cpg_map_to_hg19"}'  && \
rm ${CPGMAP_7WAY_DIR}/chr*.7orth.stitched.cpg_map_hg19
for ass in `echo hg19 panTro4 gorGor3 rheMac3 mm10 rn5 canFam3`; do
  prefix=${CPGMAP_7WAY_DIR}/${ass}_cpg_map_to_hg19
  sort -k6,6 -k7,7g -u ${prefix} -o ${prefix}.sorted_hg19
  if [ $ass != "hg19" ]; then
    sort -k1,1 -k2,2g -u ${prefix}.sorted_hg19 -o ${prefix}.sorted_${ass}
  fi
  awk 'BEGIN{OFS="\t"} {split($1,a, ".");
    split($6,b,"."); print a[2],$2,$3,$5, b[2], $7,$7+1,$10}
' < ${prefix}.sorted_${ass} > ${prefix}.cpgmap
rm ${prefix} ${prefix}.sorted_hg19 ${prefix}.sorted_${ass}
done

