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

module load miniconda3
source activate gtdbtk_2.0.0
GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/databases/GTDB-TK/release207

cd $workfolder/allgroups

while read line; do
read genome cluster_95 domain score cluster_95_final <<< $(echo $line)
grep -w $genome GMbC.rep97.cleanbins_list | awk -v cl=$cluster_95_final '{print $1"\t"cl}'
done < GMbC.dRep_cluster95_representatives.tsv | tee GMbC.representatives_list

gtdbtk classify_wf --batchfile GMbC.representatives_list --cpus $SLURM_CPUS_PER_TASK --full_tree -x fa --pplacer_cpus $SLURM_CPUS_PER_TASK --out_dir gtdbtk_r207_out --prefix GMbC


conda deactivate



source activate mashtee_env
cut -f 1 GMbC.representatives_list > GMbC.representatives_genomes
mashtree --numcpus 24 --file-of-files GMbC.representatives_genomes --outtree GMbC.mashtree.tre
