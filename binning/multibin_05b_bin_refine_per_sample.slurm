#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_05b_refine"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

source activate biscoreto_env

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
if [ -e "$workfolder/samples/${s}/${s}.refined.out" ]; then exit; fi

cd $workfolder/samples/${s}

grep $s $workfolder/subgroups/${this_local}/VAMB_out/clusters.tsv > ${s}.vamb.contigs_to_bin.tsv

awk '{print $1"\t"$2"\tvamb"}'  ${s}.vamb.contigs_to_bin.tsv > ${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tconcoct"}'  ${s}.concoct.contigs_to_bin.tsv >> ${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tmetabat2"}'  ${s}.metabat2.contigs_to_bin.tsv >> ${s}.contigs_to_bin.tsv
awk '{print $1"\t"$2"\tmaxbin2"}'  ${s}.maxbin2.contigs_to_bin.tsv >> ${s}.contigs_to_bin.tsv

Rscript /work_ifs/sukmb276/github/BiScoReTo/BiScoReTo.R -i ${s}.contigs_to_bin.tsv --hmm ${s}.hmm -o ${s}

source activate binning_env
mkdir -p cleanbins
cat ${s}.refined.contig_to_bin.out | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $datafolder/Assembly_${s}/${s}.contigs.fa /dev/stdin --output_path cleanbins
