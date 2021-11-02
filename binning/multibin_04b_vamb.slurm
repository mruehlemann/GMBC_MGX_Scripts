#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=10gb
#SBATCH --time=10-00:00
#SBATCH --output=/home/sukmb276/Isilon/Metagenomes/projects/210916_GMBC_MGX/log/%A_%a.out

###########################
####    SETUP     #######
###########################

workfolder="/work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_metagenome_data/"
source activate metagenome_env

all_local=($(cut -f 3 metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

echo $SLURMD_NODENAME
####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' $workfolder/metadata_gmbc_bn10_complete.tsv)
cd $workfolder/subgroups/${this_local}

mkdir -p binning/

###########################
####    VAMB        #######
###########################
cd /tmp
touch commands.txt
echo "module load miniconda3" > commands.txt
echo "source activate vamb_gpu_env" >> commands.txt
echo "vamb --cuda --outdir $workfolder/subgroups/${this_local}/binning/VAMB_out --fasta $workfolder/subgroups/${this_local}/${this_local}.catalogue.fna.gz --jgi $workfolder/subgroups/${this_local}/${this_local}.allsamples_depth.txt -o GMBC --minfasta 200000" >> commands.txt

ssh sukmb276@medg01 < commands.txt

cd $workfolder/subgroups/${this_local}

grep ">" binning/VAMB_out/bins/*.fna | cut -d '/' -f 4 | sed 's/:>/\ /' > binning/vamb.contig_to_bin.txt
rm -r binning/VAMB_out
