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

echo $s
mkdir -p $workfolder/samples/${s}

cat $workfolder/samples/${s}/${s}.depth.txt | awk -v s=$s '{if(NR>1){print s"_"$0}else{print }}' > ${s}.depth.txt

### metabat

mkdir -p metabat
metabat2 -i $datafolder/Assembly_${s}/$s.contigs.fa -a $workfolder/samples/${s}/${s}.depth.txt -o ${s}_metabat/${s}.metabat2 -t ${SLURM_CPUS_PER_TASK} -m 2000
for fa in $(ls ${s}_metabat/${s}.metabat2.*.fa); do mv $fa ${fa}sta; done

grep '>' ${s}_metabat/*fasta | cut -d '/' -f 2 | sed 's/:>/\ /' > $workfolder/samples/${s}/${s}.metabat2.contigs_to_bin.tsv

### maxbin2_bins
module load gcc

awk '{if(NR>1) print $1"\t"$3}' ${s}.depth.txt > ${s}.abu

mkdir -p ${s}_maxbin2_out
perl /work_ifs/sukmb276/software/MaxBin-2.2.4/run_MaxBin_prodigal.pl -contig $datafolder/Assembly_${s}/$s.contigs.fa -prodigal $workfolder/samples/${s}/${s}.prodigal.faa -out ${s}_maxbin2_out/${s}.maxbin2 -abund ${s}.abu -thread ${SLURM_CPUS_PER_TASK} -min_contig_length 2000

grep '>' ${s}_maxbin2_out/*fasta | cut -d '/' -f 2 | sed 's/:>/\ /' > $workfolder/samples/${s}/${s}.maxbin2.contigs_to_bin.tsv

#### concoct
source activate binning_env

cut_up_fasta.py $datafolder/Assembly_${s}/$s.contigs.fa -c 10000 -o 0 --merge_last -b $s.minimap.bed > ${s}.filtered.10k.fna

cp $datafolder/Assembly_${s}/$s.sorted.bam .
samtools index ${s}.sorted.bam

cat $s.minimap.bed | sed "s/${s}\_k141/k141/g" > $s.minimap.renamed.bed
concoct_coverage_table.py $s.minimap.renamed.bed $s.sorted.bam > ${s}.coverage_table.tsv
cat ${s}.coverage_table.tsv | sed "s/k141/${s}_k141/g" > ${s}.coverage_table.renamed.tsv
concoct --composition_file ${s}.filtered.10k.fna --coverage_file ${s}.coverage_table.renamed.tsv -c 1000 -r 151 -t ${SLURM_CPUS_PER_TASK} -l 2000 -s 1234 -i 500 -b ${s}

merge_cutup_clustering.py ${s}_clustering_gt2000.csv > ${s}_clustering_merged.csv

awk -F',' -v sample=$s '{if(NR>1) print sample"_concoct_bin_"$2".fasta\t"$1}'  ${s}_clustering_merged.csv > $workfolder/samples/${s}/${s}.concoct.contigs_to_bin.tsv

#done
