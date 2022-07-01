#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_07e_quant"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

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
ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/${s}.salmon.out" ]; then exit; fi

cd $TMPDIR

source activate metagenome_env

samtools bam2fq $datafolder/Assembly_${s}/$s.sorted.bam > $s.fastq

sed -i '/^@/ s/[.][12][/]/\//' $s.fastq

repair.sh in=$s.fastq out1=$s.R1.fastq out2=$s.R2.fastq outs=$s.RU.fastq interleaved=true overwrite=true

source activate binning_env

salmon quant -i $workfolder/allgroups/salmon_index -l IU -1 $s.R1.fastq -2 $s.R2.fastq --validateMappings -o out -p ${SLURM_CPUS_PER_TASK} --meta

cp out/quant.sf $workfolder/samples/${s}/${s}.salmon.out

cd $workfolder/samples/${s}/

module load R
Rscript  $workfolder/GMBC_MGX_Scripts/binning/multibin_07e_salmon_collect.R
