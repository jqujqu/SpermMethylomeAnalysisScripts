#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export ADS_DIR=~/panfs/tools/git/adssrc
export LC_ALL=C
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes

# merge sites to get regions (desert=1k)
SPECIES_LIST=(Human Chimp Gorilla Rhesus Mouse Rat Dog)
for i in `echo ${SPECIES_LIST[*]}`; do
bedtools merge -d 1000 -i ${i}_Sperm_hg19.methposterior > ${i}_Sperm_hg19.mergedbed
done

# collapse regions to get 7orth regions
${ADS_DIR}/collapsebed -c 7 -o 7orth.bed *_Sperm_hg19.mergedbed

# commonCpG sites
sort -k1,1 -k2,2g -k3,3g -m *_Sperm_hg19.methposterior | \
cut -f 1-3 | uniq -c | awk 'BEGIN{OFS="\t"} $1==7{print $2,$3,$4}' > commonCpG.bed
# all sites
sort -k1,1 -k2,2g -k3,3g -m -u *_Sperm_hg19.methposterior | \
awk '{OFS="\t"; print $1,$2,$2+1, "X", "0", "+"}' > Sites.bed


# 7way orth sites
bedtools intersect -u -a Sites.bed -b 7orth.bed > 7orth.sites.bed
awk '{print $1":"$2}' 7orth.sites.bed | sort -k1,1 -o 7orth.namesort
for i in `echo ${SPECIES_LIST[*]}`; do
  # sort sperm methposterior by name
  awk '{OFS="\t"; print $1":"$2,$5}' < ${i}_Sperm_hg19.methposterior | \
    sort -k1,1 > ${i}_Sperm_hg19.methposterior.namesort
  # fill in missing sites
  join -j1 -a1  -e "-1" -o "0 2.2" 7orth.namesort \
    ${i}_Sperm_hg19.methposterior.namesort \
    > ${i}_Sperm_hg19.methposterior.namesort.7orth
done

# make table
suffix=methposterior.namesort.${7orth}
paste Human_Sperm_hg19.${suffix}  \
Chimp_Sperm_hg19.${suffix}   \
Gorilla_Sperm_hg19.${suffix}  \
Rhesus_Sperm_hg19.${suffix}  \
Mouse_Sperm_hg19.${suffix}  \
Rat_Sperm_hg19.${suffix} \
Dog_Sperm_hg19.${suffix} | \
awk 'BEGIN{OFS="\t";} {split($1,a,":"); print a[1],a[2],$2,$4,$6,$8,$10,$12,$14}' > ${7orth}_sites.${suffix}
done

echo ${SPECIES_LIST[*]} | tr ' ' '\011' > 7orth_sites.methposterior.tmp
sort -k1,1 -k2,2g 7orth_sites.methposterior.namesort.7orth >> 7orth_sites.methposterior
# clean up
rm *namesort* *.mergedbed

# breakdown entire table into ~30 parts
mkdir ${PROJ_DIR}/hg19Multiz7orthMethylomes/breakdown
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes/breakdown
awk -v sp="${SPECIES_LIST[*]}" '
NR==2{
    part=0; chrom=$1; count=1; end=$2;
    print sp > "methposterior.table.part"part;
    print > "methposterior.table.part"part
}
NR>2 && count>500000 {
    if($2-end>1000 || $1!=chrom) {
        part+=1; count=0;
        print sp > "methposterior.table.part"part
    }
}
NR>2 {
    print > "methposterior.table.part"part;
    count+=1; end=$2; chrom=$1;
}
' < ../7orth_sites.methposterior

