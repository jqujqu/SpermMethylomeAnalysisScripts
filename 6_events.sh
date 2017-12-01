#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export HGSIZES=~/Annotation/hg19.chrom.sizes
export SCRIPT_DIR=${PROJ_DIR}/scripts
export LC_ALL=C

mkdir -p ${PROJ_DIR}/hg19Multiz7orthMethylomes/events/filtered
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes/events

cat > child_parent.txt <<- END
Euarchontoglires Boreoeutheria
Catarrhini Euarchontoglires
Homininae Catarrhini
Hominini Homininae
Human Hominini
Chimp Hominini
Gorilla Homininae
Rhesus Catarrhini
Murinae Boreoeutheria
Mouse Murinae
Rat Murinae
Dog Boreoeutheria
END

########################################################
# Call Events
########################################################
awk '{split($4,a,"_"); print > "7sp."a[1]".HYPO"}' < ../7orth_sites.seg.out

while read line ; do
chd=`echo $line|awk '{print $1}'`;
par=`echo $line|awk '{print $2}'`;
echo $chd
# birth
bedtools closest -a 7sp.${chd}.HYPO -b 7sp.${par}.HYPO -d | \
awk '$13 > 0 && $3-$2 >50 {print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g -u > 7sp.${chd}.HYPO.birth
# death
bedtools closest -b 7sp.${chd}.HYPO -a 7sp.${par}.HYPO -d | \
awk '$13 > 0{print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g -u > 7sp.${chd}.HYPO.death
# extension
bedtools closest -a 7sp.${chd}.HYPO -b 7sp.${par}.HYPO -d | \
awk '$13 == 0{print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g > tmp && \
bedtools subtract -a tmp -b 7sp.${par}.HYPO | \
sort -k1,1 -k2,2g -k3,3g -u >  7sp.${chd}.HYPO.extension && rm tmp
# contraction
bedtools closest -b 7sp.${chd}.HYPO -a 7sp.${par}.HYPO -d | \
awk '$13 == 0{print}' | cut -f 1-6 | \
sort -k1,1 -k2,2g -k3,3g -u > tmp && \
bedtools subtract -a tmp -b 7sp.${chd}.HYPO | \
sort -k1,1 -k2,2g -k3,3g -u >  7sp.${chd}.HYPO.contraction && rm tmp
done <child_parent.txt

########################################################
# Filter events (remove events adjacent to deserts)
########################################################
# deserts -- complement of the 7-way orthologous genome
${PROJ_DIR}/scripts/complement_bed.py ${HGSIZES} ../7orth.bed | \
tr ' ' '\011' < /dev/stdin > ../7orth_complement.bed

for i in `ls *HYPO.birth *HYPO.extension *HYPO.death *HYPO.contraction `; do
bedtools closest -a $i -b ../7orth_complement.bed -d | \
awk '$10>1{print}' | cut -f 1-6 > filtered/$(basename $i)
done

########################################################
# Summarize event number and size
########################################################
# (total size (Mbp), number, average_size (bp)) x
# (birth, extension, death, contraction)
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes/events/filtered
ln -s ../child_parent.txt

while read line ; do
chd=`echo $line|awk '{print $1}'`;
b=`awk '{i+=$3-$2} END{print i/1e6, NR, i/(NR+0.01)}' < 7sp.${chd}.HYPO.birth`;
e=`awk '{i+=$3-$2} END{print i/1e6, NR, i/(NR+0.01)}' < 7sp.${chd}.HYPO.extension`;
d=`awk '{i+=$3-$2} END{print i/1e6, NR, i/(NR+0.01)}' < 7sp.${chd}.HYPO.death`;
c=`awk '{i+=$3-$2} END{print i/1e6, NR, i/(NR+0.01)}' < 7sp.${chd}.HYPO.contraction`;
echo $b $e $d $c
done < child_parent.txt > tmp_summary.txt

paste child_parent.txt tmp_summary.txt | \
tr ' ' '\011' < /dev/stdin | \
cut -f 1,3-14 > events_number_size && rm tmp_summary.txt

# plot (Figure2D)
Rscript ${PROJ_DIR}/scripts/plot_event_size.R
