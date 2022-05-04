#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10Gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

source activate metagenome_env

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    SAMPLE SELECTION     ########
####################################

## selects all samples that belong to the respective subgroup (locality = column 3 of the metadata sheet)
this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples_indices=$(awk -veth=${this_local} '{if($3==eth) print (NR-1)}' /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | tr '\n' ',' | sed 's/[,]$//')

##################################
####    start slurms     ########
##################################

sbatch --array=${local_samples_indices}%10 /work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/GMBC_Scripts/binning/multibin_07b_virus_pipeline_unbinned.slurm
sbatch --array=${local_samples_indices}%10 /work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/GMBC_Scripts/binning/multibin_07c_virus_pipeline_cleanbins.slurm
