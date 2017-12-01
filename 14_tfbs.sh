#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export GALAXY_DIR=~/galaxy
export RMAP_DIR=~/rmapbs/bin
export SCRIPT_DIR=${PROJ_DIR}/scripts
export LC_ALL=C

cd ${PROJ_DIR}/TFBS/regions
EVENT_DIR=${PROJ_DIR}/hg19Multiz7orthMethylomes/events
sort -k1,1 -k2,2g -k3,3g -m ${EVENT_DIR}/7sp.Human.HYPO ${EVENT_DIR}/7sp.Chimp.HYPO \
${EVENT_DIR}/7sp.Gorilla.HYPO ${EVENT_DIR}/7sp.Rhesus.HYPO | \
${RMAP_DIR}/${RMAP_DIR}/sortbed -c /dev/stdin -o Primate_species_mergedHYPO.bed

sort -k1,1 -k2,2g -k3,3g -m ${EVENT_DIR}/7sp.Mouse.HYPO ${EVENT_DIR}/7sp.Rat.HYPO | \
${RMAP_DIR}/sortbed -c /dev/stdin -o Rodent_species_mergedHYPO.bed

# merge events on each lineage and require HYPO in extant species
for event in `echo birth extension`; do
sort -k1,1 -k2,2g -k3,3g -m ${EVENT_DIR}/filtered/7sp.Catarrhini.HYPO.${event} \
${EVENT_DIR}/filtered/7sp.Homininae.HYPO.${event} \
${EVENT_DIR}/filtered/7sp.Hominini.HYPO.${event} \
${EVENT_DIR}/filtered/7sp.Human.HYPO.${event} | \
${RMAP_DIR}/sortbed -c /dev/stdin -o Primate_lineage.${event} && \
bedtools intersect -a ${EVENT_DIR}/7sp.Human.HYPO -b Primate_lineage.${event}  | \
sort -k1,1 -k2,2g -k3,3g -u -o tmp &&  mv tmp Primate_lineage.${event}
done

for event in `echo birth extension`; do
sort -k1,1 -k2,2g -k3,3g -m ${EVENT_DIR}/filtered/7sp.Murinae.HYPO.${event} \
${EVENT_DIR}/filtered/7sp.Mouse.HYPO.${event} | \
${RMAP_DIR}/sortbed -c /dev/stdin -o Rodent_lineage.${event} && \
bedtools intersect -a ${EVENT_DIR}/7sp.Mouse.HYPO -b Rodent_lineage.${event}  | \
sort -k1,1 -k2,2g -k3,3g -u -o tmp &&  mv tmp Rodent_lineage.${event}
done

# Remove overlap between lineages
for event in `echo birth extension`; do
bedtools subtract -a Primate_lineage.${event} -b Rodent_lineage.birth > tmp && \
bedtools subtract -a tmp -b Rodent_lineage.extension | \
sort -k1,1 -k2,2g -k3,3g -u | \
${RMAP_DIR}/mapsifter -e -t Rodent_species_mergedHYPO.bed /dev/stdin | \
awk '$3-$2>100{print}' > Primate_lineage_specific.${event} && rm tmp

bedtools subtract -a Rodent_lineage.${event} -b Primate_lineage.birth > tmp && \
bedtools subtract -a tmp -b Primate_lineage.extension | \
sort -k1,1 -k2,2g -k3,3g -u | \
${RMAP_DIR}/mapsifter -e -t Primate_species_mergedHYPO.bed /dev/stdin | \
awk '$3-$2>100{print}' > Rodent_lineage_specific.${event} && rm tmp
done


# roimethstat
hg19DIR=${PROJ_DIR}/hg19Multiz7orthMethylomes
for event in `echo birth extension`; do
roimethstat -L Rodent_lineage_specific.${event} ${hg19DIR}/Mouse_Sperm_hg19.meth \
-o Rodent_lineage_specific.${event}.mouse_roi
roimethstat -L Primate_lineage_specific.${event} ${hg19DIR}/Human_Sperm_hg19.meth \
-o Primate_lineage_specific.${event}.human_roi
done

# filter: minimum 3CpG
for event in `echo birth extension`; do
awk '{split($4,a, ":"); } a[2]>=3{print}' \
< Rodent_lineage_specific.${event}.mouse_roi \
> Rodent_lineage_specific.${event}.filtered
awk '{split($4,a, ":"); } a[2]>=3{print}' \
< Primate_lineage_specific.${event}.human_roi \
> Primate_lineage_specific.${event}.filtered
done

human_hmr=${PROJ_DIR}/hg19Multiz7orthMethylomes/Human_Sperm_hg19.hmr
mouse_hmr=${PROJ_DIR}/hg19Multiz7orthMethylomes/Mouse_Sperm_hg19.hmr
TSS=~/Annotation/EnsemblGenes_hg19.TSS.bed.unique

# identify promoter HMR extensions
${RMAP_DIR}/mapsifter ${human_hmr} -t ${TSS} | \
${RMAP_DIR}/mapsifter Primate_lineage_specific.extension.filtered -t /dev/stdin -o Primate_lineage_specific.extension.filtered.at_promoter
${RMAP_DIR}/mapsifter ${mouse_hmr} -t ${TSS} | \
${RMAP_DIR}/mapsifter Rodent_lineage_specific.extension.filtered -t /dev/stdin -o Rodent_lineage_specific.extension.filtered.at_promoter
# collect the rest of gain events
sort -k1,1 -k2,2g -m Primate_lineage_specific.*.filtered | \
${RMAP_DIR}/mapsifter -e -t Primate_lineage_specific.extension.filtered.at_promoter /dev/stdin \
-o Primate_lineage_specific.gain.filtered.non_promoter
sort -k1,1 -k2,2g -m Rodent_lineage_specific.*.filtered | \
${RMAP_DIR}/mapsifter -e -t Rodent_lineage_specific.extension.filtered.at_promoter /dev/stdin \
-o Rodent_lineage_specific.gain.filtered.non_promoter


##################################################
# TFBS:
# mouse (mm9)
# https://www.encodeproject.org/metadata/type=experiment&replicates.library.biosample.donor.organism.scientific_name=Mus%20musculus&assay_term_name=ChIP-seq&target.investigated_as=transcription%20factor&files.file_type=bed%20narrowPeak/metadata.tsv
# peaks merged by target, and liftOver to hg19
# ${PROJ_DIR}/TFBS/TFBS_ENCODE_Mouse_hg19/<TF>-mouse.merged.bed.lift2hg19

# human (hg19):
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClusteredV3.bed.gz
# peaks merged by target
# ${PROJ_DIR}/TFBS/TFBS_ENCODE_Human_hg19/<TF>.bed

##################################################
mkdir -p ${PROJ_DIR}/TFBS/human_lineage_counts
cd ${PROJ_DIR}/TFBS/human_lineage_counts

REGION_DIR=${PROJ_DIR}/TFBS/regions
HumanHMR=${PROJ_DIR}/hg19Multiz7orthMethylomes/Human_Sperm_hg19.hmr
HumanTF_DIR=${PROJ_DIR}/TFBS/TFBS_ENCODE_Human_hg19

for region in `ls ${HumanHMR} \
${REGION_DIR}/Primate_lineage_specific.birth.filtered \
${REGION_DIR}/Primate_lineage_specific.extension.filtered \
${REGION_DIR}/Primate_lineage_specific.gain.filtered.non_promoter \
${REGION_DIR}/Primate_lineage_specific.extension.filtered.at_promoter`;
do REG=`basename $region `;
echo $REG;
for tf in `ls ${HumanTF_DIR}/*.bed `;
do TF=`basename $tf`;
echo $TF; NUM=`${RMAP_DIR}/mapsifter $tf -t $region | wc -l`; echo $TF $NUM >> ${REG}.TFcount ;
done;
done

# plot (Figure 6A)
Rscript tfbs_human_plot.R

##################################################
mkdir -p ${PROJ_DIR}/TFBS/mouse_lineage_counts
cd ${PROJ_DIR}/TFBS/mouse_lineage_counts

REGION_DIR=${PROJ_DIR}/TFBS/regions
MouseHMR=${PROJ_DIR}/hg19Multiz7orthMethylomes/Mouse_Sperm_hg19.hmr
MouseTF_DIR=${PROJ_DIR}/TFBS/TFBS_ENCODE_Mouse_hg19

for region in `ls ${MouseHMR} \
${REGION_DIR}/Rodent_lineage_specific.birth.filtered \
${REGION_DIR}/Rodent_lineage_specific.extension.filtered \
${REGION_DIR}/Rodent_lineage_specific.extension.filtered.at_promoter \
${REGION_DIR}/Rodent_lineage_specific.gain.filtered.non_promoter`;
do REG=`basename $region `;
echo $REG;
for tf in `ls ${MouseTF_DIR}/*.lift2hg19 `;
do TF=`basename $tf | awk '{split($0, a, "."); print $1}'`;
echo $TF; NUM=`${RMAP_DIR}/mapsifter $tf -t $region | wc -l`;
echo $TF $NUM >> ${REG}.MouseTFcount ;
done;
done

# plot (Figure 6A)
Rscript ${SCRIPT_DIR}/tfbs_mouse_plot.R

