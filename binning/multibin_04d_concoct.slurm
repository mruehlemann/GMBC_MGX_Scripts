#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
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

source activate binning_env
module load samtools/1.9

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    GROUP SELECTION     ########
####################################

## selects all samples that belong to the respective subgroup (locality = column 3 of the metadata sheet)
this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv)
cd $TMPDIR


#############################
####    concoct       #######
#############################


for s in ${local_samples}; do
#ava=$(grep $s $workfolder/samples_available.txt)
ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then continue; fi
echo $s
mkdir -p $workfolder/samples/${s}
cut_up_fasta.py $datafolder/Assembly_${s}/$s.contigs.fa -c 10000 -o 0 --merge_last -b $s.minimap.bed > ${s}.filtered.10k.fna
cp $datafolder/Assembly_${s}/$s.sorted.bam .
samtools index ${s}.sorted.bam
concoct_coverage_table.py $s.minimap.bed $s.sorted.bam > ${s}.coverage_table.tsv
concoct --composition_file ${s}.filtered.10k.fna --coverage_file ${s}.coverage_table.tsv -c 1000 -r 151 -t ${SLURM_CPUS_PER_TASK} -l 2000 -s 1234 -i 500 -b ${s}

merge_cutup_clustering.py ${s}_clustering_gt2000.csv > ${s}_clustering_merged.csv

awk -F',' -v sample=$s '{if(NR>1) print sample"_concoct_bin_"$2".fasta\t"$1}'  ${s}_clustering_merged.csv > $workfolder/samples/${s}/${s}.concoct.contigs_to_bin.tsv

done
