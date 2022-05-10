#!/bin/bash
#SBATCH -c 6
#SBATCH --mem=60Gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_06b_checkm"

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
all_samples=($(ls $workfolder/samples))
#all_samples=($(ls -d $datafolder/Assembly_* | awk -F'_' '{print $NF}'))
s=${all_samples[$SLURM_ARRAY_TASK_ID]}

this_local=$(grep -w $s /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | cut -f 3 )

#############################
####    MAXBIN DEFAULT #######
#############################
cd $TMPDIR

source activate phamb-1.0.2_env
PHAMB_DIR="/work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources"

ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/virus_detection/unbinned_contigs/${s}.checkv.out" ]; then exit; fi

cat $workfolder/samples/${s}/cleanbins/*.fa > allbinned.fa

python3 $PHAMB_DIR/DeepVirFinder/dvf.py -i allbinned.fa -o DVF -l 2000

checkv contamination allbinned.fa checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv completeness allbinned.fa checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv complete_genomes allbinned.fa checkv_out
checkv quality_summary allbinned.fa checkv_out

mkdir -p $workfolder/samples/${s}/virus_detection/binned_contigs/phamb_output

cp checkv_out/quality_summary.tsv $workfolder/samples/${s}/virus_detection/binned_contigs/${s}.checkv.binned.out
cat checkv_out/*.fna > $workfolder/samples/${s}/virus_detection/binned_contigs/${s}.checkv.candidates.fna
