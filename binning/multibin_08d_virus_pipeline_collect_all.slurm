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

source activate phamb-1.0.2_env
PHAMB_DIR="/work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources"

#ava=$(ls -d $datafolder/Assembly_${s})
#if [ "$ava" == "" ]; then exit; fi
#if [ -e "$workfolder/samples/${s}/virus_detection/unbinned_contigs/${s}.checkv_final.out" ]; then exit; fi

cat $workfolder/samples/*/virus_detection/*.viral_final.fna > GMbC.viral_final.fna

mkdir -p $workfolder/allgroups/viral_catalog/blast_dbs
cp GMbC.viral_final.fna  $workfolder/allgroups/viral_catalog/

### from checkv website: https://bitbucket.org/berkeleylab/checkv/src/master/
#First, create a blast+ database:
makeblastdb -in GMbC.viral_final.fna -dbtype nucl -out  $workfolder/allgroups/viral_catalog/blast_dbs/GMbC.viral.db

cat GMbC.viral_final.fna  /work_ifs/ikmb_repository/shared/microbiome/databases/MGV_v1.0_2021_07_08/mgv_votu_representatives.fna > GMbC_plus_MGV.viral.fna
makeblastdb -in GMbC_plus_MGV.viral.fna -dbtype nucl -out  $workfolder/allgroups/viral_catalog/blast_dbs/GMbC_plus_MGV.viral.db


# #Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
# blastn -query all.viral_final.fna -db all.viral.db -outfmt '6 std qlen slen' -max_target_seqs 10000 -out viral.blast.tsv -num_threads ${SLURM_CPUS_PER_TASK}
#
# #Note: using the -perc_identity flag will speed up the search at the cost of sensitivity:
# #blastn -query <my_seqs.fna> -db <my_db> -outfmt '6 std qlen slen' -max_target_seqs 10000 -perc_identity 90 -o <my_blast.tsv> -num_threads 32
#
# #Next, calculate pairwise ANI by combining local alignments between sequence pairs:
# python $PHAMB_DIR/checkv/scripts/anicalc.py -i viral.blast.tsv -o all.ani.tsv
#
# #Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
# python $PHAMB_DIR/checkv/scripts/aniclust.py --fna all.viral_final.fna --ani all.ani.tsv --out all.viral_cluster.tsv --min_ani 95 --min_tcov 85 --min_qcov 0
#
