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
source activate binning_env

cd $workfolder/subgroups/$this_local/
mkdir -p genome_links

for i in $local_samples; do ls $workfolder/samples/$i/cleanbins/*.fa; done | awk -F '/' '{print $0"\t"$NF}' | sed 's/[.]fa$//' >  $workfolder/subgroups/$this_local/$this_local.cleanbins_list

head -n 1 $workfolder/samples/$(echo $local_samples | cut -d ' ' -f 1)/*.refined.out > ${this_local}.refined.out
for i in $local_samples; do cat $workfolder/samples/$i/$i.refined.out | grep -v "uniqueSCGs.gtdb_rel207_bac120"; done >> ${this_local}.refined.out

cd $TMPDIR
cut -f 1 $workfolder/subgroups/$this_local/$this_local.cleanbins_list > genomes_list
dRep compare -g genomes_list -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -sa 0.95 -d .
cp data_tables/Cdb.csv ${this_local}.dRep_fastANI_95.csv
python $workfolder/GMBC_MGX_Scripts/binning/multibin_07a_drep_fastANI99.py
cp data_tables/Cdb_97.csv ${this_local}.dRep_fastANI_97.csv
cp data_tables/Cdb_99.csv ${this_local}.dRep_fastANI_99.csv
cp ${this_local}.dRep_fastANI_95.csv ${this_local}.dRep_fastANI_97.csv ${this_local}.dRep_fastANI_99.csv $workfolder/subgroups/$this_local
conda deactivate

cd $workfolder/subgroups/$this_local/
module load R
Rscript  $workfolder/GMBC_MGX_Scripts/binning/multibin_07a_drep_get_representative.R
