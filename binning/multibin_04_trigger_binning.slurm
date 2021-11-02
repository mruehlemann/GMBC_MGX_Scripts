#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10gb
#SBATCH --time=10-00:00
#SBATCH --output=/home/sukmb276/Isilon/Metagenomes/projects/210916_GMBC_MGX/log/%A_%a.out

###########################
####    SETUP     #######
###########################

workfolder="/work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_metagenome_data/"
#source activate metagenome_env

all_local=($(cut -f 3 metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

echo $SLURMD_NODENAME
####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' $workfolder/metadata_gmbc_bn10_complete.tsv)
cd $workfolder


##################################
####    start slurms     ########
##################################

sbatch --array=${SLURM_ARRAY_TASK_ID} GMBC_Scripts/binning/multibin_04a_metabat.slurm
sbatch --array=${SLURM_ARRAY_TASK_ID} GMBC_Scripts/binning/multibin_04c_maxbin_ar.slurm
sbatch --array=${SLURM_ARRAY_TASK_ID} GMBC_Scripts/binning/multibin_04c_maxbin_bac.slurm
sbatch --array=${SLURM_ARRAY_TASK_ID} GMBC_Scripts/binning/multibin_04c_maxbin_def.slurm
#sbatch --array=${SLURM_ARRAY_TASK_ID} GMBC_Scripts/binning/multibin_04b_vamb.slurm
