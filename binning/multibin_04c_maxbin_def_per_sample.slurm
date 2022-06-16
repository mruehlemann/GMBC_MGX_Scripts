#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_04c_maxbin"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

source activate metagenome_env
module load gcc

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    SAMPLE SELECTION     ########
####################################

## selects all samples that belong to the respective subgroup (locality = column 3 of the metadata sheet)
#this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
#mkdir -p $workfolder/subgroups/${this_local}
#local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv
all_samples=($(cut -f 1 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv))
#all_samples=($(ls -d $datafolder/Assembly_* | awk -F'_' '{print $NF}'))
s=${all_samples[$SLURM_ARRAY_TASK_ID]}
cd $TMPDIR


#############################
####    MAXBIN DEFAULT #######
#############################
ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/${s}.maxbin2.contigs_to_bin.tsv" ]; then exit; fi

awk '{if(NR>1) print $1"\t"$3}' $workfolder/samples/${s}/${s}.depth.txt > ${s}.abu

mkdir -p ${s}_maxbin2_out
perl /work_ifs/sukmb276/software/MaxBin-2.2.4/run_MaxBin_prodigal.pl -contig $datafolder/Assembly_${s}/$s.contigs.fa -prodigal $workfolder/samples/${s}/${s}.prodigal.faa -out ${s}_maxbin2_out/${s}.maxbin2 -abund ${s}.abu -thread ${SLURM_CPUS_PER_TASK} -min_contig_length 2000

grep '>' ${s}_maxbin2_out/*fasta | cut -d '/' -f 2 | sed 's/:>/\ /' > $workfolder/samples/${s}/${s}.maxbin2.contigs_to_bin.tsv
rm -r binning/maxbin_multi_out/
