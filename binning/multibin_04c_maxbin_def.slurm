#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
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

mkdir -p binning/maxbin_multi_out/
perl /work_ifs/sukmb276/software/MaxBin-2.2.4/run_MaxBin_prodigal.pl -contig $this_local.catalogue.fna.gz -prodigal ${this_local}.prodigal.faa -out binning/maxbin_multi_out/${this_local}.DEFAULTMULTI -abund_list ${this_local}.abulist -thread ${SLURM_CPUS_PER_TASK} -min_contig_length 2000

grep '>' binning/maxbin_multi_out/*fasta | cut -d '/' -f 3 | sed 's/:>/\ /' > binning/maxbin_default.contig_to_bin.txt
rm -r binning/maxbin_multi_out/
