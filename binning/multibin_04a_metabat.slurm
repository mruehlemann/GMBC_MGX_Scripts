#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=240gb
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


###########################
####    METABAT     #######
###########################

mkdir -p binning/metabat2_multi_out/
metabat2 -i ${this_local}.catalogue.fna.gz -a  ${this_local}.allsamples_depth.txt -o binning/metabat2_multi_out/${this_local}.METABAT2MULTI -t ${SLURM_CPUS_PER_TASK} -m 2000
for fa in $(ls binning/metabat2_multi_out/${this_local}.METABAT2MULTI.*.fa); do mv $fa ${fa}sta; done

grep '>' binning/metabat2_multi_out/*fasta | cut -d '/' -f 3 | sed 's/:>/\ /' > binning/metabat.contig_to_bin.txt
rm -r binning/metabat2_multi_out
