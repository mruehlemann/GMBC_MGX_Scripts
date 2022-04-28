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

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv)
cd $TMPDIR

####################################
####    CATALOG CREATION    ########
####################################

for s in ${local_samples}; do
#ava=$(grep $s $workfolder/samples_available.txt)
ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then continue; fi
echo $s
cat $datafolder/Assembly_${s}/$s.contigs.fa | cut -f 1 -d ' ' | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >= 2000 {print ">"$0}' > ${s}.sample.fasta
done

cat *.sample.fasta > ${this_local}.catalogue.fna

minimap2 -I100G -d ${this_local}.catalogue.mmi ${this_local}.catalogue.fna # make index

gzip ${this_local}.catalogue.fna

cp ${this_local}.catalogue.mmi ${this_local}.catalogue.fna.gz  $workfolder/subgroups/${this_local}
