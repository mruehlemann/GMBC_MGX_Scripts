#!/bin/bash
#SBATCH -c 2
#SBATCH --mem=20Gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_06b_checkm"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    SAMPLE SELECTION     ########
####################################

## selects all samples that belong to the respective subgroup (locality = column 3 of the metadata sheet)
#this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
#mkdir -p $workfolder/subgroups/${this_local}
#local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv)
all_samples=($(cut -f 1 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv))
#all_samples=($(ls -d $datafolder/Assembly_* | awk -F'_' '{print $NF}'))
#s=${all_samples[$SLURM_ARRAY_TASK_ID]}

#this_local=$(grep -w $s /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | cut -f 3 )

#############################
####    MAXBIN DEFAULT #######
#############################
cd $TMPDIR

source activate vContact2_env
PHAMB_DIR="/work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources"

cat $workfolder/samples/*/virus_detection/unbinned_contigs/*.viral_final.fna > all.viral_final.fna
s="all"
mkdir -p ${s}_tmp_workfolder/
cat ${s}.viral_final.fna | parallel -j ${SLURM_CPUS_PER_TASK} --block 999k --recstart '>' --pipe prodigal -p meta -a ${s}_tmp_workfolder/${s}.{#}.faa -d ${s}_tmp_workfolder/${s}.{#}.ffn -o tmpfile
cat ${s}_tmp_workfolder/$s.*.faa > ${s}.faa

echo "protein_id,contig_id,keywords" > $s.g2g.csv
grep '>' ${s}.faa | tr -d '>' | awk '{print $1","$1","$1}' | sed -r 's/\_[0-9]+$//' | sed -r 's/\_viral\_[0-9]+\_[0-9]+//' | awk -F ',' '{print $2","$3","$1}' >> $s.g2g.csv

vcontact2 --raw-proteins ${s}.faa --rel-mode 'Diamond' --proteins-fp $s.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/sukmb276/.conda/envs/vContact2_env/bin/cluster_one-1.0.jar --output-dir ${s}_vcontact2_out
