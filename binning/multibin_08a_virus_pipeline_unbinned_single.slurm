#!/bin/bash
#SBATCH -c 6
#SBATCH --mem=60Gb
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
s=${all_samples[$SLURM_ARRAY_TASK_ID]}

this_local=$(grep -w $s /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | cut -f 3 )

#############################
####   PHAMB         #######
#############################
#
# This pipleine follows the PhAMB pipeline, however a bit smarter, as contigs binned by the refinement procedure are
# not included in the processing => PhAMB works based on bins and would classify these as bacterial anyway
# The clean bins are processed separately with an approach not relying on PhAMB
# In addition, not only the VAMB input, but also CONCOCT and METABAT2 output are used as sources for potential viral bins

cd $TMPDIR

source activate phamb-1.0.2_env
PHAMB_DIR="/work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources"

mkdir -p $workfolder/samples/${s}/virus_detection/single_contigs

ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/virus_detection/single_contigs/${s}.checkv.out" ]; then exit; fi

mkdir phamb_discovery

### Identification of short and binned contigs for exclusion
grep '>' $workfolder/samples/${s}/cleanbins/*fa | cut -d '>' -f 2 | cut -d ' ' -f 1 > phamb_discovery/bact_contigs.txt
cat $datafolder/Assembly_${s}/${s}.contigs.fa | cut -d ' ' -f 1 | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) < 10000 {print $1"\n"}' > phamb_discovery/short_contigs.txt
cat phamb_discovery/bact_contigs.txt phamb_discovery/short_contigs.txt | cut -d ' ' -f 1 > phamb_discovery/bad_contigs.txt

/work_ifs/sukmb276/software/bin/faSomeRecords $datafolder/Assembly_${s}/${s}.contigs.fa -exclude phamb_discovery/bad_contigs.txt /dev/stdout | cut -d ' ' -f 1 > phamb_discovery/contigs.fna

### Gene calls from unbinned contigs are retained for annotation
grep '>' phamb_discovery/contigs.fna | tr -d '>' | cut -f 1 -d ' ' > phamb_discovery/good_contigs.txt
#cat phamb_discovery/bad_contigs.txt | cut -f 1 | awk '{for(i=1;i<=1000;++i) print $1"_"i}' | /work_ifs/sukmb276/software/bin/faSomeRecords $workfolder/samples/${s}/${s}.prodigal.faa -exclude /dev/stdin phamb_discovery/proteins.faa
cat phamb_discovery/good_contigs.txt | cut -f 1 | awk '{for(i=1;i<=5000;++i) print $1"_"i}' | /work_ifs/sukmb276/software/bin/faSomeRecords $workfolder/samples/${s}/${s}.prodigal.faa /dev/stdin phamb_discovery/proteins.faa


hmmsearch -Z 1 --cpu ${SLURM_CPUS_PER_TASK} --noali --tblout imgvr.out /work_ifs/ikmb_repository/shared/microbiome/databases/MGV/viral_detection_pipeline/input/imgvr.hmm phamb_discovery/proteins.faa
hmmsearch -Z 1 --cut_tc --cpu ${SLURM_CPUS_PER_TASK} --noali --tblout pfam.out /work_ifs/ikmb_repository/shared/microbiome/databases/MGV/viral_detection_pipeline/input/pfam.hmm phamb_discovery/proteins.faa

ln -s /work_ifs/ikmb_repository/shared/microbiome/databases/MGV/viral_detection_pipeline/input/ input

python /work_ifs/ikmb_repository/shared/microbiome/databases/MGV/viral_detection_pipeline/count_hmm_hits.py phamb_discovery/contigs.fna phamb_discovery/proteins.faa imgvr.out pfam.out > hmm_hits.tsv
python /work_ifs/ikmb_repository/shared/microbiome/databases/MGV/viral_detection_pipeline/strand_switch.py phamb_discovery/contigs.fna phamb_discovery/proteins.faa > strand_switch.tsv

conda activate r_microbiome_env
Rscript /work_ifs/ikmb_repository/shared/microbiome/databases/MGV/viral_detection_pipeline/virfinder.R phamb_discovery/contigs.fna virfinder.tsv

source activate phamb-1.0.2_env

python /work_ifs/ikmb_repository/shared/microbiome/databases/MGV/viral_detection_pipeline/master_table.py hmm_hits.tsv virfinder.tsv strand_switch.tsv > master_table.tsv

python /work_ifs/ikmb_repository/shared/microbiome/databases/MGV/viral_detection_pipeline/viral_classify_short.py \
--features master_table.tsv \
--out_base ${s}.viral

cp ${s}.viral.tsv $workfolder/samples/${s}/virus_detection/single_contigs/

cat ${s}.viral.tsv | cut -f 1 | /work_ifs/sukmb276/software/bin/faSomeRecords phamb_discovery/contigs.fna /dev/stdin $s.candidates.fna

checkv contamination ${s}.candidates.fna checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv completeness ${s}.candidates.fna checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv complete_genomes ${s}.candidates.fna checkv_out
checkv quality_summary ${s}.candidates.fna checkv_out

cp checkv_out/quality_summary.tsv $workfolder/samples/${s}/virus_detection/single_contigs/${s}.checkv.out
cat checkv_out/*.fna > $workfolder/samples/${s}/virus_detection/single_contigs/${s}.checkv.candidates.fna
