#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=120Gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_02b_map"

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

this_local=$(grep -w $s /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | cut -f 3 )

###########################
####    MAPPING     #######
###########################

mkdir -p bams

#for s in ${local_samples}; do
#ava=$(grep $s $workfolder/samples_available.txt)
ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then continue; fi
if [ -e "$workfolder/subgroups/${this_local}/${this_local}.${s}.depth.txt" ]; then exit; fi
echo $s
minimap2 -t ${SLURM_CPUS_PER_TASK} -N 50 -ax sr  $workfolder/subgroups/${this_local}/${this_local}.catalogue.mmi $datafolder_reads/${s}_P_1.fastq.gz $datafolder_reads/${s}_P_2.fastq.gz > $s.minimap.sam
samtools view -F 3584 -b --threads ${SLURM_CPUS_PER_TASK} $s.minimap.sam | samtools sort > bams/$s.minimap.bam
rm $s.minimap.sam
#done

jgi_summarize_bam_contig_depths --outputDepth $workfolder/subgroups/${this_local}/${this_local}.${s}.depth.txt bams/*
