#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=100gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_01_catalog"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
source activate metagenome_env

###########################
####    all subgroups     #######
###########################

allcluster=($(cat $workfolder/allgroups/GMbC.dRep_cluster95_representatives.tsv | cut -f 5))
cluster=${allcluster[$SLURM_ARRAY_TASK_ID]}

tax="s__Prevotella copri_A"
tax2=$(echo $tax | tr ' ' '_')
mkdir -p cophylo_analysis

cd $TMPDIR
mkdir $tax2
cd $tax2
parallel -j 8 prodigal -i $workfolder/samples/{= s:\_cleanbin\_[0-9]+[.]fa:: =}/cleanbins/{} -a $(echo "{}a") ::: $(grep -w "$tax" $workfolder/allgroups/GMbC_final_tax.tsv | cut -f 1 | xargs -I {} grep {} $workfolder/allgroups/GMbC.all_genomes.tsv | cut -f 1)

export PATH=$PATH:/work_ifs/sukmb276/software/bin

python /work_ifs/sukmb276/github/CoreCruncher/corecruncher_master.py -in . -out cc_out -freq 70 -align mafft
