#!/bin/bash
#SBATCH --reservation=sukmb412_13
#SBATCH --gres=gpu:1
#SBACTH --nodelist=medg01
#SBATCH -c 16
#SBATCH --mem=500gb
#SBATCH --time=3-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_04b_vamb"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

#source activate metagenome_env

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))


###########################
####    VAMB        #######
###########################

source activate metagenome_env
cd $TMPDIR

while read line; do
echo $line
read path genome <<< $(echo $line)
awk -v genome=$genome '/^>/ {print ">"genome"_contig_" ++i; next}{print}' $path > ${genome}.fna
done < $workfolder/allgroups/GMbC.representatives_list

cat GMbC_cluster95_*.fna >  GMbC.catalogue.fna

minimap2 -I100G -d GMbC.catalogue.mmi GMbC.catalogue.fna # make index

conda deactivate
conda activate binning_env

salmon index -t GMbC.catalogue.fna -i GMbC_salmon_index

gzip GMbC.catalogue.fna

cp GMbC.catalogue.mmi GMbC.catalogue.fna.gz  $workfolder/subgroups/allgroups
cp -r GMbC_salmon_index $workfolder/subgroups/allgroups


#
#
# source activate mashtee_env
# cut -f 1 GMbC.representatives_list > GMbC.representatives_genomes
# mashtree --numcpus 24 --file-of-files GMbC.representatives_genomes --outtree GMbC.mashtree.tre
