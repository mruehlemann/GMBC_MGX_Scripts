#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_05x_dastool"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

module load R

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    SAMPLE SELECTION     ########
####################################

## selects all samples that belong to the respective subgroup (locality = column 3 of the metadata sheet)
#this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
#mkdir -p $workfolder/subgroups/${this_local}
#local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv)
all_samples=($(cut -f 1 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv))
#all_samples=($(ls -d $datafolder/Assembly_* | awk -F'_' '{print $NF}'))
s=${all_samples[$SLURM_ARRAY_TASK_ID]}

this_local=$(grep -w $s /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | cut -f 3 )


#############################
####    MAXBIN DEFAULT #######
#############################

source activate dastool_env
module load ruby

export PATH=$PATH:/work_ifs/sukmb276/software/DAS_Tool/

cd /work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/samples/${s}
mkdir dastool_input
awk '{print $2"\t"$1}' ${s}.maxbin2.contigs_to_bin.tsv > dastool_input/maxbin2.out
awk '{print $2"\t"$1}' ${s}.metabat2.contigs_to_bin.tsv > dastool_input/metabat2.out
awk '{print $2"\t"$1}' ${s}.concoct.contigs_to_bin.tsv > dastool_input/concoct.out
awk '{print $2"\t"$1}' ${s}.vamb.contigs_to_bin.tsv > dastool_input/vamb.out

mkdir dastool_out
DAS_Tool -i dastool_input/maxbin2.out,dastool_input/metabat2.out,dastool_input/concoct.out,dastool_input/vamb.out -c $datafolder/Assembly_${s}/${s}.contigs.fa -o dastool_out/${s}.dastool -p ${s}.prodigal.faa --write_bin_evals --labels=maxbin,metabat,concoct,vamb -t 8 --write_bins

rm -r dastool_input
