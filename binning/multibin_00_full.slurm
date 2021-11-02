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
if [ "$ava" == "" ]; then next; fi
echo $s
zcat $workfolder/assemblies/$s/$s.contigs.fasta.gz  | sed 's/Spades_//' | sed 's/_contig//' | awk -F'_' '/^>/ {printf $1"GMBC%.10i\n", $2; next} {print}' > ${s}.sample.fasta
done

concatenate.py --keepnames -m 2000 ${this_local}.catalogue.fna.gz *.sample.fasta

minimap2 -I100G -d ${this_local}.catalogue.mmi ${this_local}.catalogue.fna.gz; # make index

cp ${this_local}.catalogue.mmi ${this_local}.catalogue.fna.gz  $workfolder/subgroups/${this_local}


###########################
####    MAPPING     #######
###########################

mkdir -p bams

for s in ${local_samples}; do
ava=$(grep $s $workfolder/samples_available.txt)
if [ "$ava" == "" ]; then next; fi
echo $s
minimap2 -t ${SLURM_CPUS_PER_TASK} -N 50 -ax sr ${this_local}.catalogue.mmi $datafolder/${s}_P_1.fastq.gz $datafolder/${s}_P_2.fastq.gz > $s.minimap.sam
samtools view -F 3584 -b --threads 20 $s.minimap.sam | samtools sort > bams/$s.minimap.bam
rm $s.minimap.sam
done

jgi_summarize_bam_contig_depths --outputDepth $workfolder/subgroups/${this_local}/${this_local}.allsamples_depth.txt bams/


######################################################
######### CDS CALLING AND ANNOTATION     ############
######################################################

mkdir -p tmp_workfolder
zcat ${this_local}.catalogue.fna.gz | parallel -j ${SLURM_CPUS_PER_TASK} --block 999k --recstart '>' --pipe prodigal -p meta -a tmp_workfolder/${this_local}.{#}.faa -d tmp_workfolder/${this_local}.{#}.fna -o tmpfile
cat tmp_workfolder/$this_local.*.faa > $workfolder/subgroups/${this_local}/${this_local}.prodigal.faa
cat tmp_workfolder/$this_local.*.fna > $workfolder/subgroups/${this_local}/${this_local}.prodigal.fna

source activate gtdbtk_1.5.0

GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/shared/microbiome/databases/GTDB/release202/

cat $workfolder/subgroups/${this_local}/${this_local}.prodigal.faa > TEST_protein.faa
python ~/Isilon/software/annot.py
grep -v "^#" TEST_ANNOT_OUT/TEST/TEST_pfam.tsv | awk '{print $1" - "$6" "$6" "$13" "$12" 0 "$13" "$12" 1 1 0 0 1 1 1 1 "}' | sort -k 3,3 > ${this_local}.annot.pfam
grep -v "^#" TEST_ANNOT_OUT/TEST/TEST_tigrfam.tsv > ${this_local}.annot.tigrfam
cat ${this_local}.annot.tigrfam  ${zshort}.annot.pfam | sed -r 's/[ ]+/ /g' > ${this_local}.annot

rm -r TEST_ANNOT_OUT tmp_workfolder tmpfile

for i in $(cat ~/Isilon/software/binning_helpers/markers_bac120_clean.txt); do grep $i ${this_local}.annot; done > ${this_local}.annot.bac120
for i in $(cat ~/Isilon/software/binning_helpers/markers_ar122_clean.txt); do grep $i ${this_local}.annot; done > ${this_local}.annot.ar122

cp ${this_local}.annot ${this_local}.annot.bac120 ${this_local}.annot.ar122
