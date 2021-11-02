#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/home/sukmb276/Isilon/Metagenomes/projects/210916_GMBC_MGX/log/%A_%a.out

###########################
####    SETUP     #######
###########################

workfolder="/work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_metagenome_data/"
source activate metagenome_env
module load gcc

all_local=($(cut -f 3 metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

echo $SLURMD_NODENAME
####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' $workfolder/metadata_gmbc_bn10_complete.tsv)
cd $workfolder/subgroups/${this_local}


#############################
####    MAXBIN DEFAULT #######
#############################

module load R

cd binning

cat *.contig_to_bin.txt > contig_to_bin.txt

Rscript /work_ifs/sukmb276/github/DAS_Tool_ripoff_new_param.Rscript $this_local contig_to_bin.txt ../${this_local}.annot

mkdir -p ../cleanbins

for i in $(cut -f 1 ${this_local}.refined.contig_to_bin.out | sort | uniq | grep -v binnew); do
grep $i ${this_local}.refined.contig_to_bin.out | cut -f 2 | /work_ifs/sukmb276/software/bin/faSomeRecords ../${this_local}.catalogue.fna.gz /dev/stdin ../cleanbins/$i
done
