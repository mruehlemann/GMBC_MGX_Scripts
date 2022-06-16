#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_04a_metabat"

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
#this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
#mkdir -p $workfolder/subgroups/${this_local}
#local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv)
all_samples=($(cut -f 1 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv))
#all_samples=($(ls -d $datafolder/Assembly_* | awk -F'_' '{print $NF}'))
s=${all_samples[$SLURM_ARRAY_TASK_ID]}
cd $TMPDIR

###########################
####    METABAT     #######
###########################

#for s in ${local_samples}; do
#ava=$(grep $s $workfolder/samples_available.txt)
ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/${s}.metabat2.contigs_to_bin.tsv" ]; then exit; fi
echo $s
mkdir -p $workfolder/samples/${s}
mkdir -p metabat
metabat2 -i $datafolder/Assembly_${s}/$s.contigs.fa -a $workfolder/samples/${s}/${s}.depth.txt -o ${s}_metabat/${s}.metabat2 -t ${SLURM_CPUS_PER_TASK} -m 2000
for fa in $(ls ${s}_metabat/${s}.metabat2.*.fa); do mv $fa ${fa}sta; done

grep '>' ${s}_metabat/*fasta | cut -d '/' -f 2 | sed 's/:>/\ /' > $workfolder/samples/${s}/${s}.metabat2.contigs_to_bin.tsv
