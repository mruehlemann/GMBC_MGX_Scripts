#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40Gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_02b_map"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

#workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
workfolder="/work_ifs/sukmb276/Metagenomes/projects/220512_GMbC_NovaSeq/temporary_workdir/"

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
cd $TMPDIR

this_local=$(grep -w $s /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | cut -f 3 )

###########################
####    MAPPING     #######
###########################

module load miniconda3
source activate metaphlan2_env
export PERL5LIB=

#for s in ${local_samples}; do
#ava=$(grep $s $workfolder/samples_available.txt)
ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then continue; fi
mkdir -p $workfolder/samples/$s/
reads=$(ls $datafolder_reads/${s}_P_1.fastq.gz)
if [ "$reads" == "" ]; then

module load samtools/1.9

samtools bam2fq $datafolder/Assembly_${s}/$s.sorted.bam > $s.fastq
sed -i '/^@/ s/[.][12][/]/\//' $s.fastq
repair.sh in=$s.fastq out1=$s.R1.fastq out2=$s.R2.fastq outs=$s.RU.fastq interleaved=true overwrite=true
cat $s.R1.fastq $s.R2.fastq | metaphlan2.py --input_type fastq --nproc $SLURM_CPUS_PER_TASK -x mpa_v20_m200 --bowtie2db /work_ifs/ikmb_repository/databases/Metaphlan/2.0/ > $workfolder/samples/$s/$s.metaphlan2.out

else
zcat $datafolder_reads/${s}_P_1.fastq.gz $datafolder_reads/${s}_P_2.fastq.gz | metaphlan2.py --input_type fastq --nproc $SLURM_CPUS_PER_TASK -x mpa_v20_m200 --bowtie2db /work_ifs/ikmb_repository/databases/Metaphlan/2.0/ > $workfolder/samples/$s/$s.metaphlan2.out

fi
