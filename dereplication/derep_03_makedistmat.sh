#!/bin/bash

# Based on Bacsort script

workfolder="/work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_metagenome_data/"
source activate metagenome_env

~/Isilon/github/Bacsort/scripts/pairwise_identities_to_distance_matrix.py --max_dist 0.2 genome_collection/fastani_output > genome_collection/fastani.phylip
cat genome_collection/fastani.phylip | awk -F'/' '/^\// {print $NF; next}; {print}' | sed 's/[.]fasta//' > genome_collection/fastani_final.phylip
cat genome_collection/mash.phylip | awk -F'/' '/^\// {print $NF; next}; {print}' | sed 's/[.]fasta//' > genome_collection/mash_final.phylip
cd genome_collection
~/Isilon/github/Bacsort/scripts/combine_distance_matrices.py fastani_final.phylip mash_final.phylip > final.phylip
