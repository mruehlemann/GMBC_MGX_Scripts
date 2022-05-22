#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=1-00:00
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

####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' $workfolder/metadata_gmbc_bn10_with_novaseq.tsv)

###########################
####    VAMB        #######
###########################

module load miniconda3
source activate gtdbtk_2.1.0
GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/databases/GTDB-TK/release207_v2

cd $TMPDIR
#cd $workfolder/allgroups

cut -f 1 $workfolder/subgroups/${this_local}/${this_local}.dRep_cluster97_representatives.tsv | xargs -I {} grep -w {} $workfolder/subgroups/${this_local}/${this_local}.cleanbins_list > batch

gtdbtk classify_wf --batchfile batch --cpus $SLURM_CPUS_PER_TASK -x fa --pplacer_cpus $SLURM_CPUS_PER_TASK --out_dir gtdbtk_r207v2_out --prefix GMbC_${this_local}

conda deactivate

cp gtdbtk_r207v2_out/GMbC_${this_local}.*.summary.tsv $workfolder/subgroups/${this_local}
