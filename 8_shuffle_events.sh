#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export HGSIZES=~/Annotation/hg19.chrom.sizes
export SCRIPT_DIR=${PROJ_DIR}/scripts
mkdir -p ${PROJ_DIR}/hg19Multiz7orthMethylomes/events/filtered/shuffle_events
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes/events/filtered/shuffle_events
export LC_ALL=C

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

while read line ; do
echo $line
child=`echo $line | awk '{print $1}'`
parent=`echo $line | awk '{print $2}'`

sort -k1,1 -k2,2g -m ../7sp.${child}.HYPO.extension \
../7sp.${child}.HYPO.birth -o tmp_gain

bedtools intersect -a ../../7sp.${child}.HYPO \
-b ../../7sp.${parent}.HYPO | \
sort -k1,1 -k2,2g -k3,3g -u > tmp_conserved

sort -k1,1 -k2,2g -m ../7sp.${child}.HYPO.contraction \
../7sp.${child}.HYPO.death -o tmp_loss

cp ../../7sp.${parent}.HYPO tmp_anchypo
bedtools subtract -a ../../../7orth.bed -b tmp_anchypo |\
sort -k1,1 -k2,2g -k3,3g -u > tmp_anchyper

bedtools closest -a tmp_gain -b tmp_conserved -d | \
cut -f 13 > branch_${child}_gain_dist_conserved.txt
bedtools closest -a tmp_loss -b tmp_conserved -d | \
cut -f 13 > branch_${child}_loss_dist_conserved.txt

for i in `seq 1 100`; do echo $i
bedtools shuffle -noOverlapping -i tmp_gain -g ${HGSIZES} -incl tmp_anchyper | \
sort -k1,1 -k2,2g -k3,3g -u | awk 'NF==6{print}' > tmp_gain_shuffled
bedtools closest -a tmp_gain_shuffled -b tmp_conserved -d | \
cut -f 13 >> branch_${child}_gain_dist_conserved_shuffled.txt
done

for i in `seq 1 100`; do echo $i
bedtools shuffle -noOverlapping -i tmp_loss -g ${HGSIZES} -incl tmp_anchypo | \
sort -k1,1 -k2,2g -k3,3g -u | awk 'NF==6{print}' >  tmp_loss_shuffled
bedtools closest -a tmp_loss_shuffled -b tmp_conserved -d | \
cut -f 13 >> branch_${child}_loss_dist_conserved_shuffled.txt
done
rm tmp*
done < child_parent.txt

# plot (Figure S6)
Rscript shuffle_events_plot.R
