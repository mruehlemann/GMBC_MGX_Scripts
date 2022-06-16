#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/home/sukmb276/Isilon/Metagenomes/projects/210916_GMBC_MGX/log/%A_%a.out

# Based on Bacsort script

workfolder="/work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_metagenome_data/"
source activate metagenome_env

mash sketch -l genome_collection/mq.paths.txt -p 24 -o genome_collection/reference -s 100000

mash dist -p 24 genome_collection/reference.msh genome_collection/reference.msh -t > genome_collection/distances.tab

# Convert matrix to PHYLIP format
tail -n +2 genome_collection/distances.tab > genome_collection/distances.tab.temp  # remove first line
wc -l genome_collection/distances.tab.temp | awk '"[0-9]+ errors" {sum += $1}; END {print sum}' > genome_collection/distances.ndist  # get number for sample/line count
cat genome_collection/distances.ndist genome_collection/distances.tab.temp > genome_collection/mash.phylip
rm genome_collection/distances.ndist genome_collection/distances.tab.temp
printf "\n"
