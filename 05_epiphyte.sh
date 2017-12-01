#!/bin/bash

export PROJ_DIR=~/sperm_methylome_evolution
export EPIPHYTE_DIR=~/Epiphyte/bin
export LC_ALL=C
SPECIES_LIST=(Human Chimp Gorilla Rhesus Mouse Rat Dog)

########################################################
# Estimation
########################################################
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes/breakdown
# starting point of branch lengths
echo "(((((Human:0.007,Chimp:0.007)HC:0.003,Gorilla:0.01)HCG:0.02,Rhesus:0.03)Primate:0.06,(Mouse:0.02,Rat:0.02)Rodent:0.07)PR:0.01,Dog:0.1)Root;" > 7sp_binary.nwk

# estimate phylo-epigenetic tree with missing data
for part in `seq 0 32`; do
${EPIPHYTE_DIR}/epiphy-est \
-h 1000 -f -v -o part${part}.params \
7sp_binary.nwk methposterior.table.part${part} &> log_part${part}
done

# parameter medians from 33 genomic segments:
#-------------------
#$ cat median_params
#tree (((((Human:0.438242,Chimp:0.4834855)Hominini:0.174121,Gorilla:0.5339045)Homininae:0.522948,Rhesus:0.6575265)Catarrhini:0.5655645,(Mouse:0.80988,Rat:0.8654615)Murinae:0.97688)Euarchontoglires:0.0224654,Dog:0.8801645)Boreoeutheria:0;
#pi0     0.0779114
#rate0   0.922704
#g0[0]   0.99574
#g1[0]   0.998957
#g0[1]   0.9869015
#g1[1]   0.8929415
#-------------------

# Impute methylation probabilities at missing sites;
# Estimate branch lengths and mutation rates
# using independent-site model
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes
for i in `seq 0 6`; do
SPECIES=${SPECIES_LIST[i]}; echo $SPECIES
awk -v i=$i ' BEGIN{OFS="\t"}
    NR>1 {print $1,$2,$2+1, $(i+3)}
    ' < 7orth_sites.methposterior > tmp.0
bedtools closest -a tmp.0 \
    -b ${PROJ_DIR}/hg19Multiz7orthMethylomes/${SPECIES}_Sperm_hg19.methposterior \
    -fu -D "ref" > tmp.u
bedtools closest -a tmp.0 \
    -b ${PROJ_DIR}/hg19Multiz7orthMethylomes/${SPECIES}_Sperm_hg19.methposterior \
    -fd -D "ref" > tmp.d
paste tmp.u tmp.d | \
    awk 'BEGIN{OFS="\t"}
        $4!=-1 {print $1,$2, $4}
        $4==-1 && ($11<= -1000 || $22>=1000) {print $1,$2, -1}
        $4==-1 && $11> -1000 && $22<1000 {
            print $1, $2, ($9*(-$11)+$20*($22))/($22-$11)
        }
' > tmp.${SPECIES}.interpolate && rm tmp.0 tmp.u tmp.d
done

echo ${SPECIES_LIST[*]} > 7orth_sites.methposterior.interpolate
paste tmp.Human.interpolate tmp.Chimp.interpolate \
    tmp.Gorilla.interpolate tmp.Rhesus.interpolate \
    tmp.Mouse.interpolate tmp.Rat.interpolate \
    tmp.Dog.interpolate | \
    cut -f 1-3,6,9,12,15,18,21 >> 7orth_sites.methposterior.interpolate \
    && rm tmp.*.interpolate

${EPIPHYTE_DIR}/epiphy-est -f -b 50 -k 20 -h 100  -i 300 \
    -d 0 -v -o indep.params breakdown/7sp_binary.nwk  \
    7orth_sites.methposterior.interpolate  2> indeplog

#-------------------
#$ cat indep.params
#tree    (((((Human:0.0386967,Chimp:0.0264527)Hominini:0.009452,Gorilla:0.0407904)HCG:0.0208771,Rhesus:0.0726092)Primate:0.0288665,(Mouse:0.0737205,Rat:0.113226)Rodent:0.0572602)Supraprimate:0.0196637,Dog:0.0580591)Root:0;
#pi0     0.136774
#rate0   0.509019
#-------------------


########################################################
# inference
########################################################
cd ${PROJ_DIR}/hg19Multiz7orthMethylomes

# parameters used for ancestral methylome reconstruction
cat > model.params <<-END
tree    (((((Human:0.0386967,Chimp:0.0264527)Hominini:0.009452,Gorilla:0.0407904)Homininae:0.0208771,Rhesus:0.0726092)Catarrhini:0.0288665,(Mouse:0.0737205,Rat:0.113226)Murinae:0.0572602)Euarchontoglires:0.0196637,Dog:0.0580591)Boreoeutheria:0;
pi0    0.0779114
rate0    0.509019
g0[0]    0.99574
g1[0]    0.998957
g0[1]    0.95
g1[1]    0.95
END

# posterior methylation probabilities
${EPIPHYTE_DIR}/epiphy-post \
-h 10000 -b 10000 -v -o 7orth_sites.post.out \
model.params 7orth_sites.methposterior  &> postlog

# segmentation
${EPIPHYTE_DIR}/epiphy-seg -v -o 7orth_sites.seg.out \
model.params 7orth_sites.post.out


