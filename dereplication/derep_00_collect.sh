#!/bin/bash

workfolder="/work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_metagenome_data/"
source activate metagenome_env

mkdir -p genome_collection
cat subgroups/*/*.checkm.out | grep -v "Bin ID" > genome_collection/all.checkm.out
cat genome_collection/all.checkm.out | awk -F '\t' '{if($12 >= 50 && $13 < 10) print}' > genome_collection/mq.checkm.out
for i in $(cut -f 1 genome_collection/mq.checkm.out); do
pop=$(echo $i | sed -r 's/_cleanbin_[0-9]+//');
echo $workfolder/subgroups/$pop/cleanbins/${i}.fasta
done | tee genome_collection/mq.paths.txt
