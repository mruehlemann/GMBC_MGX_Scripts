#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_07c_gtdb"

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
source activate gtdbtk_2.1.0
GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/databases/GTDB-TK/release207_v2

cd $TMPDIR
#cd $workfolder/allgroups
awk -vbatch=$SLURM_ARRAY_TASK_ID '{if(NR % 10 == batch) print}' $workfolder/allgroups/GMbC.representatives_list > batch

gtdbtk classify_wf --batchfile batch --cpus $SLURM_CPUS_PER_TASK -x fa --pplacer_cpus $SLURM_CPUS_PER_TASK --out_dir gtdbtk_r207v2_out --prefix GMbC_batch_$SLURM_ARRAY_TASK_ID

conda deactivate

cp gtdbtk_r207v2_out/GMbC_batch_${SLURM_ARRAY_TASK_ID}.*.summary.tsv $workfolder/allgroups
