#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80Gb
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
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' $workfolder/metadata_gmbc_bn10_with_novaseq.tsv)

#############################
####    MAXBIN DEFAULT #######
#############################
cd $TMPDIR

source activate vContact2_env
PHAMB_DIR="/work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources"

for s in ${local_samples}; do
cat $workfolder/samples/${s}/virus_detection/unbinned_contigs/*.viral_final.fna > ${s}.viral_final.fna
done

cat *.viral_final.fna > ${this_local}.viral.fna

mkdir -p ${this_local}_tmp_workfolder/
cat ${this_local}.viral.fna | parallel -j ${SLURM_CPUS_PER_TASK} --block 999k --recstart '>' --pipe prodigal -p meta -a ${this_local}_tmp_workfolder/${this_local}.{#}.faa -d ${this_local}_tmp_workfolder/${this_local}.{#}.ffn -o tmpfile
cat ${this_local}_tmp_workfolder/${this_local}.*.faa > ${this_local}.faa

#cp ${this_local}.faa ${this_local}.plus_millard.faa
#zcat /work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources/vcontact2/millard_phages/vcontact_format_30May2020.faa.gz  >> ${this_local}.plus_millard.faa
cat ${this_local}.faa /work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources/vcontact2/inphared_db/1Aug2022_vConTACT2_proteins.faa  > ${this_local}.plus_inphared.faa


echo "protein_id,contig_id,keywords" > ${this_local}.g2g.csv
grep '>' ${this_local}.faa | tr -d '>' | awk '{print $1","$1","$1}' | sed -r 's/\_[0-9]+$//' | sed -r 's/\_viral\_[0-9]+\_[0-9]+//' | awk -F ',' '{print $2","$3","$1}' >> ${this_local}.g2g.csv
#cp ${this_local}.g2g.csv ${this_local}.plus_millard.g2g.csv
#zcat /work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources/vcontact2/millard_phages/gene-to-genome.30May2020.csv.gz >> ${this_local}.plus_millard.g2g.csv
cat ${this_local}.g2g.csv /work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources/vcontact2/inphared_db/1Aug2022_vConTACT2_gene_to_genome.csv > ${this_local}.plus_inphared.g2g.csv


#vcontact2 --raw-proteins ${this_local}.plus_millard.faa --rel-mode 'Diamond' --proteins-fp $this_local.g2g.csv --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/sukmb276/.conda/envs/vContact2_env/bin/cluster_one-1.0.jar --output-dir ${this_local}_vcontact2_out
vcontact2 --raw-proteins ${this_local}.plus_inphared.faa --rel-mode 'Diamond' --proteins-fp $this_local.plus_inphared.g2g.csv --db 'None' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin /home/sukmb276/.conda/envs/vContact2_env/bin/cluster_one-1.0.jar --output-dir ${this_local}_vcontact2_out

mkdir -p $workfolder/subgroups/${this_local}/vcontact2_inphared_out
cp ${this_local}_vcontact2_out/viral_cluster_overview.csv ${this_local}_vcontact2_out/c1.ntw ${this_local}_vcontact2_out/genome_by_genome_overview.csv $workfolder/subgroups/${this_local}/vcontact2_inphared_out

## keep based on clustering taxonomy
## next step: cluster unassigned candidates across all subgroups, without database
