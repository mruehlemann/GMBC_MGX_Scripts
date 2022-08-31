#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_04b_vamb"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

#source activate metagenome_env

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' $workfolder/metadata_gmbc_bn10_with_novaseq.tsv)

###########################
####    VAMB        #######
###########################

module load miniconda3

cd $workfolder/subgroups/$this_local/

source activate phamb-1.0.2_env
cd /work_ifs/ikmb_repository/shared/microbiome/databases/MGV/crispr_spacers

~/Isilon/software/bin/parallel -j ${SLURM_CPUS_PER_TASK} --colsep "\t" python identify_crispr.py -i {1} -o $TMPDIR/{2} :::: $workfolder/subgroups/$this_local/$this_local.cleanbins_list

cd $TMPDIR
for i in $(ls -d *); do cat $i/*.spacers | grep -E "[[:space:]][ATGC]+$" | awk -vi=$i '{if(length($NF)>=20 && length($NF)<=80) print i"\t"$NF}' | grep -v spacer_seq | sort | uniq ; done >  $workfolder/subgroups/$this_local/$this_local.crispr_spacers.tsv
