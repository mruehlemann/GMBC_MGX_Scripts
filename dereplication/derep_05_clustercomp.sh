#!/bin/bash

# Based on Bacsort script

workfolder="/work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_metagenome_data/"
source activate metagenome_env

cd genome_collection

for thisbin in $(cut -f 4 fastANI_cluster_AvLinkage_refined.tsv | sort | uniq | grep -v "^$"); do

genomes=$(awk -vbin=$thisbin '{if($4==bin) print $1}' fastANI_cluster_AvLinkage_refined.tsv | tr -d '\"' | sed 's/[.]L/-L/g')

for g in $genomes; do
sty=$(echo $g)
#cat $targ | sed 's/[.:-]/_/g' > clusters/FGB/FGB_${thisbin}/fasta/${sty}_$g
echo ${thisbin} $sty
done | sort | uniq -c | awk '{print $2, $3, $1}'  |  tr -d '\"'

done | tee FGB_cluster_orig.txt

###

for thisbin in $(cut -f 3 fastANI_cluster_AvLinkage_refined.tsv | sort | uniq); do

genomes=$(awk -vbin=$thisbin '{if($3==bin) print $1}' fastANI_cluster_AvLinkage_refined.tsv | tr -d '\"'  | sed 's/[.]L/-L/g')

for g in $genomes; do
targ=$(cat ../01_Candidates/list_MAGs.txt | grep "$g.fasta$")
sty=$(grep $(echo $g | cut -d '-' -f 1) ../01_Candidates/smp_orig.txt | cut -f 1)
if [ "$sty" == "" ]; then sty="Ref"; fi
#cat $targ | sed 's/[-:.]/_/g' > clusters/SGB/SGB_${thisbin}/fasta/${sty}_$g
echo ${thisbin} $sty
done | sort | uniq -c | awk '{print $2, $3, $1}'  |  tr -d '\"'

done | tee GGB_cluster_orig.txt


###

for thisbin in $(cut -f 2 fastANI_cluster_AvLinkage_refined.tsv | sort | uniq); do

genomes=$(awk -vbin=$thisbin '{if($2==bin) print $1}' fastANI_cluster_AvLinkage_refined.tsv | tr -d '\"'  | sed 's/[.]L/-L/g')

for g in $genomes; do
targ=$(cat ../01_Candidates/list_MAGs.txt | grep "$g.fasta$")
sty=$(grep $(echo $g | cut -d '-' -f 1) ../01_Candidates/smp_orig.txt | cut -f 1)
if [ "$sty" == "" ]; then sty="Ref"; fi
#cat $targ | sed 's/[-:.]/_/g' > clusters/SGB/SGB_${thisbin}/fasta/${sty}_$g
echo ${thisbin} $sty
done | sort | uniq -c | awk '{print $2, $3, $1}'  |  tr -d '\"'

done | tee SGB_cluster_orig.txt
