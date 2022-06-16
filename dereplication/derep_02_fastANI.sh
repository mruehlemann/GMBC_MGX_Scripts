#!/bin/bash

# Based on Bacsort script

workfolder="/work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_metagenome_data/"
source activate metagenome_env

group_count=100

printf "\n"
echo "Find pairwise distances with FastANI"
echo "------------------------------------------------"
mkdir -p clusters
cd clusters
cat ../genome_collection/mq.paths.txt > cluster_list

mkdir -p ../genome_collection/tmp

total_count=$( wc -l < cluster_list )
count_per_file=$( perl -w -e "use POSIX; print ceil($total_count/$group_count), qq{\n}" )
cat cluster_list | shuf > split.tmp
split -a 4 -dl $count_per_file split.tmp cluster_list_
rm split.tmp


for q in $(ls cluster_list_* | tail -n 92)*; do
    q_num=$(echo $q | sed 's/cluster_list_//')
    for r in cluster_list_*; do
        r_num=$(echo $r | sed 's/cluster_list_//')
        sbatch --nodes=1 --job-name=fastANI_"$q_num"_"$r_num" --ntasks=1 --cpus-per-task=1 --mem=4096 --time=0-24:0:00 --wrap "fastANI --rl "$r" --ql "$q" -o ../genome_collection/tmp/fastani_output_"$q_num"_"$r_num" &> ../genome_collection/tmp/fastani_stdout_"$q_num"_"$r_num
    done
    if [ "$(echo $q_num | awk '{print $1%2}')" -eq "1" ]; then echo "$q_num"; echo "waiting 180 seconds"; sleep 180; fi
done

cd $workfolder

cat genome_collection/tmp/fastani_output_* > genome_collection/fastani_output
rm -r 
