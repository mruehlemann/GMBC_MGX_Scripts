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
#cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

#source activate metagenome_env

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    GROUP SELECTION     ########
####################################

mkdir -p $workfolder/allgroups

###########################
####    VAMB        #######
###########################

module load miniconda3
source activate binning_env

head -n 1 $workfolder/samples/$(ls $workfolder/samples | head -n 1)/*.refined.out > $workfolder/allgroups/all.rep97.refined.out

for this_local in ${all_local[@]};do

cd $workfolder/subgroups/$this_local/

cut -f 1 ${this_local}.dRep_cluster97_representatives.tsv  | xargs -I {} grep {} ${this_local}.cleanbins_list >>  $workfolder/allgroups/GMbC.rep97.cleanbins_list
cut -f 1 ${this_local}.dRep_cluster97_representatives.tsv  | xargs -I {} grep {} ${this_local}.refined.out  >>  $workfolder/allgroups/GMbC.rep97.refined.out

done


cd $TMPDIR
cut -f 1 $workfolder/allgroups/GMbC.rep97.cleanbins_list | sed '/^$/d' > genomes_list
dRep compare -g genomes_list -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -sa 0.95 -d .
cp data_tables/Cdb.csv $workfolder/allgroups/GMbC.rep97.dRep_95.csv
python $workfolder/GMBC_MGX_Scripts/binning/multibin_07a_drep_fastANI99.py
cp data_tables/Cdb_97.csv $workfolder/allgroups/GMbC.rep97.dRep_97.csv
cp data_tables/Cdb_99.csv $workfolder/allgroups/GMbC.rep97.dRep_99.csv
conda deactivate

cd /$workfolder/allgroups/

module load R
Rscript  $workfolder/GMBC_MGX_Scripts/binning/multibin_07b_drep_all_get_representative.R


source activate binning_env
