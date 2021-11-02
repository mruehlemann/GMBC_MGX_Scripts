#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/home/sukmb276/Isilon/Metagenomes/projects/210916_GMBC_MGX/log/%A_%a.out

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_metagenome_data/"
source activate metagenome_env

ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' $workfolder/metadata_gmbc_bn10_complete.tsv)
cd $TMPDIR

####################################
####    CATALOG CREATION    ########
####################################

for s in ${local_samples}; do
ava=$(grep $s $workfolder/samples_available.txt)
if [ "$ava" == "" ]; then continue; fi
echo $s
zcat $workfolder/assemblies/$s/$s.contigs.fasta.gz  | sed 's/Spades_//' | sed 's/_contig//' | awk -F'_' '/^>/ {printf $1"GMBC%.10i\n", $2; next} {print}' > ${s}.sample.fasta
done

concatenate.py --keepnames -m 2000 ${this_local}.catalogue.fna.gz *.sample.fasta

minimap2 -I100G -d ${this_local}.catalogue.mmi ${this_local}.catalogue.fna.gz; # make index

cp ${this_local}.catalogue.mmi ${this_local}.catalogue.fna.gz  $workfolder/subgroups/${this_local}
