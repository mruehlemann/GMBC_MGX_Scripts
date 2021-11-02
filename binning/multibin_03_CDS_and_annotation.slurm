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
source activate gtdbtk_1.5.0

all_local=($(cut -f 3 metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

echo $SLURMD_NODENAME
####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' $workfolder/metadata_gmbc_bn10_complete.tsv)
cd $TMPDIR


######################################################
######### CDS CALLING AND ANNOTATION     ############
######################################################

mkdir -p tmp_workfolder
zcat $workfolder/subgroups/${this_local}/${this_local}.catalogue.fna.gz | ~/Isilon/software/bin/parallel -j ${SLURM_CPUS_PER_TASK} --block 999k --recstart '>' --pipe prodigal -p meta -a tmp_workfolder/${this_local}.{#}.faa -d tmp_workfolder/${this_local}.{#}.fna -o tmpfile
cat tmp_workfolder/$this_local.*.faa > $workfolder/subgroups/${this_local}/${this_local}.prodigal.faa
cat tmp_workfolder/$this_local.*.fna > $workfolder/subgroups/${this_local}/${this_local}.prodigal.fna

GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/shared/microbiome/databases/GTDB/release202/

mkdir splitfa
awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%100000==0){file=sprintf("splitfa/myseq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < $workfolder/subgroups/${this_local}/${this_local}.prodigal.faa

mkdir tmp
for i in $(ls splitfa); do
  echo $i
  cp splitfa/$i TEST_protein.faa
  python ~/Isilon/software/annot.py
  grep -v "^#" TEST_ANNOT_OUT/TEST/TEST_pfam.tsv | awk '{print $1" - "$6" "$6" "$13" "$12" 0 "$13" "$12" 1 1 0 0 1 1 1 1 "}' | sort -k 3,3 > tmp/${i}.annot.pfam
  grep -v "^#" TEST_ANNOT_OUT/TEST/TEST_tigrfam.tsv > tmp/${i}.annot.tigrfam
  rm -r TEST_ANNOT_OUT
done

cat tmp/*.annot.pfam > ${this_local}.annot.pfam
cat tmp/*.annot.tigrfam > ${this_local}.annot.tigrfam

cat ${this_local}.annot.tigrfam ${this_local}.annot.pfam | sed -r 's/[ ]+/ /g' > ${this_local}.annot

for i in $(cat ~/Isilon/software/binning_helpers/markers_bac120_clean.txt); do grep $i ${this_local}.annot; done > ${this_local}.annot.bac120
for i in $(cat ~/Isilon/software/binning_helpers/markers_ar122_clean.txt); do grep $i ${this_local}.annot; done > ${this_local}.annot.ar122

cp ${this_local}.annot ${this_local}.annot.bac120 ${this_local}.annot.ar122 $workfolder/subgroups/${this_local}/

rm -r tmp_workfolder tmpfile tmp
