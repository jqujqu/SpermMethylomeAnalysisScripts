#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export ADSSRC_DIR=~/adssrc/
export RMAP_DIR=~/rmapbs/bin
export SCRIPT_DIR=${PROJ_DIR}/scripts
export LC_ALL=C
mkdir -p ${PROJ_DIR}/GC_profile/Dog_Mouse/canFam3_profile
mkdir -p ${PROJ_DIR}/GC_profile/Dog_Mouse/mm10_profile
mkdir -p ${PROJ_DIR}/GC_profile/Human_Mouse/hg19_profile
mkdir -p ${PROJ_DIR}/GC_profile/Human_Mouse/mm10_profile

############################################################
# file locations of: TSS, genomic sequence, pair-wise aligned regions, CpG sites
HG_TSS=~/Annotation/EnsemblGenes_hg19.TSS.bed.unique
HG_CHROM_DIR=~/hg19_chroms/  ## directory containing fasta files
HG_MM_NET=~/Annotation/hg19.mm10.net.bed  ## created by ${SCRIPT_DIR}/axtNet_to_bed.sh
HG_CPG=~/Annotation/hg19_CpG.bed  ## CpG dinucleotide locations

MM_TSS=~/Annotation/EnsemblGenes_mm10.TSS.bed.unique
#directory containing fasta files
MM_CHROM=~/mm10_chroms/
MM_HG_NET=~/Annotation/mm10.hg19.net.bed
MM_CAN_NET=~/Annotation/mm10.canFam3.net.bed
MM_CPG=~/Annotation/mm10_CpG.bed

CAN_TSS=~/Annotation/EnsemblGenes_canFam3.TSS.bed.unique
#directory containing fasta files
CAN_CHROM=~/canFam3_chroms/
CAN_MM_NET=~/Annotation/canFam3.mm10.net.bed
CAN_CPG=~/Annotation/canFam3_CpG.bed

############################################################
#### OPT1
# hg19 human vs mouse
cd ${PROJ_DIR}/GC_profile/Human_Mouse/hg19_profile
MEhmr=${PROJ_DIR}/hg19Multiz100wayMethylomes/Mouse_ESC_hg19.hmr
HEhmr=${PROJ_DIR}/hg19Multiz100wayMethylomes/Human_ESC_hg19.hmr
MShmr=${PROJ_DIR}/hg19Multiz100wayMethylomes/Mouse_Sperm_hg19.hmr
HShmr=${PROJ_DIR}/hg19Multiz100wayMethylomes/Human_Sperm_hg19.hmr
TSS=${HG_TSS}
chromdir=${HG_CHROM_DIR}
Net=${HG_MM_NET}
cpg=${HG_CPG}
bedoverlap $MEhmr $HEhmr | \
${RMAP_DIR}/mapsifter $TSS -t /dev/stdin -o tmp && \
${RMAP_DIR}/mapsifter $HEhmr -t tmp -o HEhmr_TSS


#### OPT2
# mm10 human vs mouse
cd ${PROJ_DIR}/GC_profile/Human_Mouse/mm10_profile
MEhmr=${PROJ_DIR}/mm10Multiz60wayMethylomes/Mouse_ESC_mm10.hmr
HEhmr=${PROJ_DIR}/mm10Multiz60wayMethylomes/Human_ESC_mm10.hmr
MShmr=${PROJ_DIR}/mm10Multiz60wayMethylomes/Mouse_Sperm_mm10.hmr
HShmr=${PROJ_DIR}/mm10Multiz60wayMethylomes/Human_Sperm_mm10.hmr
TSS=${MM_TSS}
chromdir=${MM_CHROM}
Net=${MM_HG_NET}
cpg=${MM_CpG}
bedoverlap $MEhmr $HEhmr | \
${RMAP_DIR}/mapsifter $TSS -t /dev/stdin -o tmp && \
${RMAP_DIR}/mapsifter $HEhmr -t tmp -o HEhmr_TSS

#### OPT3
#### canFam3 dog vs mouse
cd ${PROJ_DIR}/GC_profile/Dog_Mouse/canFam3_profile

MEhmr=${PROJ_DIR}/OriginalMethylomes/Dog_canFam3/OtherSpecies/Mouse_ESC_canFam3.hmr
HEhmr=${PROJ_DIR}/OriginalMethylomes/Dog_canFam3/OtherSpecies/Mouse_ESC_canFam3.hmr
MShmr=${PROJ_DIR}/OriginalMethylomes/Dog_canFam3/OtherSpecies/Mouse_Sperm_canFam3.hmr
HShmr=${PROJ_DIR}/OriginalMethylomes/Dog_canFam3/Dog_Sperm.hmr
TSS=${CAN_TSS}
chromdir=${CAN_CHROM}
Net=${CAN_MM_NET}
cpg=${CAN_CPG}
bedoverlap $MEhmr $HEhmr | \
${RMAP_DIR}/mapsifter $TSS -t /dev/stdin -o tmp && \
${RMAP_DIR}/mapsifter $HEhmr -t tmp -o HEhmr_TSS

#### OPT4
#### mm10 dog vs mouse
cd ${PROJ_DIR}/GC_profile/Dog_Mouse/mm10_profile

MEhmr=${PROJ_DIR}/mm10Multiz60wayMethylomes/Mouse_ESC_mm10.hmr
HEhmr=${PROJ_DIR}/mm10Multiz60wayMethylomes/Mouse_ESC_mm10.hmr
MShmr=${PROJ_DIR}/mm10Multiz60wayMethylomes/Mouse_Sperm_mm10.hmr
HShmr=${PROJ_DIR}/mm10Multiz60wayMethylomes/Dog_Sperm_mm10.hmr
TSS=${MM_TSS}
chromdir=${MM_CHROM}
Net=${MM_CAN_NET}
cpg=${MM_CpG}

bedoverlap $MEhmr $HEhmr | \
${RMAP_DIR}/mapsifter $TSS -t /dev/stdin -o tmp && \
${RMAP_DIR}/mapsifter $HEhmr -t tmp -o HEhmr_TSS

############################################################
#### Common script
bedtools closest -a $MShmr -b $HEhmr -d | \
awk '$13==0{print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g -u -o tmp && \
bedtools subtract -a tmp -b $HEhmr | \
sort -k1,1 -k2,2g -k3,3g -u | \
awk '$3-$2>300{print}' > MS-HE.bed
bedtools closest -a $HShmr -b $MShmr -d | \
awk '$13==0{print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g -u -o tmp && \
bedtools subtract -a tmp -b $MShmr | \
sort -k1,1 -k2,2g -k3,3g -u | \
awk '$3-$2>300{print}' > HS-MS.bed

awk 'BEGIN{OFS="\t"} $2>500{print $1,$2-500,$2,$4,$5,$6}
{print $1,$3,$3+500,$4,$5,$6; }' < $HShmr | \
sort -k1,1 -k2,2g -k3,3g  > pad

bedtools closest -a HS-MS.bed -b pad -d | \
awk '$13==1{print}' | cut -f 1-6 | \
sortbed /dev/stdin -o tmp && mv tmp HS-MS.bed

bedtools closest -a MS-HE.bed -b HEhmr_TSS -d | \
awk '$13==1{print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g -u -o MS-HE.bed.filtered
bedtools closest -a HS-MS.bed -b MS-HE.bed.filtered -d | \
awk '$13==1{print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g -u -o HS-MS.bed.filtered

bedtools closest -a MS-HE.bed.filtered -b HS-MS.bed.filtered -d | \
awk '$13==1{print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g -u -o tmp && mv tmp MS-HE.bed.filtered
bedtools closest -a HEhmr_TSS -b MS-HE.bed.filtered -d | \
awk '$13==1{print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g -u -o HEhmr_filtered

bedtools closest -a MS-HE.bed.filtered -b HEhmr_filtered -d | \
awk 'BEGIN{OFS="\t"} $2==$9{print $1,$2,$3,$4,$5,"-"}
$3==$8{print  $1,$2,$3,$4,$5,"+"}' | \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u > MS-HE.bed.filtered.stranded

bedtools closest -a HS-MS.bed.filtered -b MS-HE.bed.filtered | \
awk 'BEGIN{OFS="\t"} $2==$9{print $1,$2,$3,$4,$5,"-"}
$3==$8{print  $1,$2,$3,$4,$5,"+"}' | \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u > HS-MS.bed.filtered.stranded

bedtools closest -a HEhmr_filtered -b MS-HE.bed.filtered -D ref | \
awk 'BEGIN{OFS="\t"} $13==1{print $1,$2,$3,$4,$5,"-" }
$13==-1{print $1,$2,$3,$4,$5,"+" }' | \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u -o HEhmr_filtered.stranded

awk 'BEGIN{OFS="\t"} $6=="+"{print $1,$2-500,$3,"ref_HSBD", $2, $6;}
$6=="-"{print $1,$2,$3+500, "ref_HSBD",$3, $6}' < HS-MS.bed.filtered.stranded | \
sort -k1,1 -k2,2g -k3,3g -k6,6 > I-II_ref_HSBD.bed

bedtools closest -a HS-MS.bed.filtered.stranded -b MS-HE.bed.filtered.stranded -d |
awk 'BEGIN{OFS="\t"} $6=="+"{print $1,$2,$9,"ref_MSBD", $3, $6}
$6=="-"{print $1,$8,$3,"ref_MSBD", $2,$6}' | \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u -o II-III_ref_MSBD.bed

bedtools closest -a MS-HE.bed.filtered.stranded -b HEhmr_filtered -d | \
awk 'BEGIN{OFS="\t"} $6=="+"{print $1,$2,$9,"ref_HEBD",$3,$6}
$6=="-"{print $1,$8,$3,"ref_HEBD", $2, $6}' | \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u -o III-IV_ref_HEBD.bed


# split II, III in half, so bases are not counted twice
region=I-II_ref_HSBD.bed
awk 'BEGIN{OFS="\t"}
$2<$5 && $3>$5 && $6=="+"{start=$2;end=int(($3+$5)/2);
if (start<$5-500) start=$5-500;
if(end>$5+500) end=$5+500;
print $1,start,end,$4,$5,$6}
$2<$5 && $3>$5 && $6=="-"{start=int(($2+$5)/2);end=$3;
if (start<$5-500) start=$5-500;
if(end>$5+500) end=$5+500;
print $1,start,end,$4,$5,$6}
' < $region > ${region}.max500bp

region=II-III_ref_MSBD.bed
awk 'BEGIN{OFS="\t"}
$2<$5 && $3>$5 {start=int(($2+$5)/2);end=int(($3+$5)/2);
if (start<$5-500) start=$5-500;
if(end>$5+500) end=$5+500;
print $1,start,end,$4,$5,$6}
' < $region > ${region}.max500bp

region=III-IV_ref_HEBD.bed
awk 'BEGIN{OFS="\t"}
$2<$5 && $3>$5 && $6=="-"{start=$2;end=int(($3+$5)/2);
if (start<$5-500) start=$5-500;
if(end>$5+500) end=$5+500;
print $1,start,end,$4,$5,$6}
$2<$5 && $3>$5 && $6=="+"{start=int(($2+$5)/2);end=$3;
if (start<$5-500) start=$5-500;
if(end>$5+500) end=$5+500;
print $1,start,end,$4,$5,$6}
' < $region > ${region}.max500bp

# make sure I doesn't overlap with other HMRs
bedtools closest  -b II-III_ref_MSBD.bed  -a $MShmr -d | \
awk '$13>0{print }' |cut -f 1-6| \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u -o  MShmr_other
bedtools closest  -b II-III_ref_MSBD.bed  -a $HShmr -d | \
awk '$13>0{print }' |cut -f 1-6| \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u -o HShmr_other
bedtools closest  -b II-III_ref_MSBD.bed  -a $MEhmr -d | \
awk '$13>0{print }' |cut -f 1-6| \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u -o  MEhmr_other
bedtools closest  -b II-III_ref_MSBD.bed  -a $HEhmr -d | \
awk '$13>0{print }' |cut -f 1-6| \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u -o HEhmr_other
sortbed -c MShmr_other HShmr_other MEhmr_other HEhmr_other -o Shmr_other
bedtools subtract -a I-II_ref_HSBD.bed.max500bp -b Shmr_other | \
sort -k1,1 -k2,2g -k3,3g -k6,6 -u -o tmp && \
mv tmp I-II_ref_HSBD.bed.max500bp

export halfwinsize=1
for region in `ls I-II_ref_HSBD.bed II-III_ref_MSBD.bed III-IV_ref_HEBD.bed`; do
echo $region
bedtools intersect -a ${region}.max500bp -b ${Net} > ${region}.max500bp.aligned
${SCRIPT_DIR}/GCprofile.py -b ${region}.max500bp.aligned \
-w ${halfwinsize} -s 1 -f ${chromdir} \
-o ${region}.max500bp.aligned.GCprofile.${halfwinsize}bp
done

############################################################
#### plot (Figure 4C,  S7C)
cd ${PROJ_DIR}/GC_profile
Rscript ${SCRIPT_DIR}/plot_gcprofile.R
