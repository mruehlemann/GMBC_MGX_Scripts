#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40Gb
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
####    MAXBIN DEFAULT #######
#############################
cd $TMPDIR

source activate phamb-1.0.2_env
PHAMB_DIR="/work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources"

ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/virus_detection/${s}.checkv_final.out" ]; then exit; fi

#cat $workfolder/samples/${s}/virus_detection/unbinned_contigs/$s.checkv.candidates.fna $workfolder/samples/${s}/virus_detection/single_contigs/$s.checkv.candidates.fna $workfolder/samples/${s}/virus_detection/binned_contigs/$s.checkv.proviruses.fna | awk '{print $1}' | sed -E 's/([._][0-9]+)([_]1$)/\1/' > ${s}.cnd.fna
cat $workfolder/samples/${s}/virus_detection/unbinned_contigs_new/$s.checkv.candidates.fna $workfolder/samples/${s}/virus_detection/binned_contigs/$s.checkv.proviruses.fna | awk '{print $1}' | sed -E 's/([._][0-9]+)([_]1$)/\1/' > ${s}.cnd.fna
### get seqs with at least low viral confidence
#cat $workfolder/samples/${s}/virus_detection/unbinned_contigs/$s.checkv.out $workfolder/samples/${s}/virus_detection/single_contigs/$s.checkv.out $workfolder/samples/${s}/virus_detection/binned_contigs/$s.checkv.proviruses.out | egrep "quality|Complete"  | awk '{if($10 >= 50 && $1!="contig_id") print}' > $s.checkv_cand_on.out
cat $workfolder/samples/${s}/virus_detection/unbinned_contigs_new/$s.checkv.out $workfolder/samples/${s}/virus_detection/binned_contigs/$s.checkv.proviruses.out | egrep "quality|Complete"  | grep -v concoct | awk '{if($10 >= 50 && $1!="contig_id") print}' > $s.checkv_cand_on.out
cat $s.checkv_cand_on.out | cut -f 1 | /work_ifs/sukmb276/software/bin/faSomeRecords ${s}.cnd.fna /dev/stdin $s.viral.fna

### from checkv website: https://bitbucket.org/berkeleylab/checkv/src/master/
#First, create a blast+ database:
makeblastdb -in $s.viral.fna -dbtype nucl -out $s.viral.db

#Next, use megablast from blast+ package to perform all-vs-all blastn of sequences:
blastn -query $s.viral.fna -db $s.viral.db -outfmt '6 std qlen slen' -max_target_seqs 10000 -out viral.blast.tsv -num_threads ${SLURM_CPUS_PER_TASK}

#Note: using the -perc_identity flag will speed up the search at the cost of sensitivity:
#blastn -query <my_seqs.fna> -db <my_db> -outfmt '6 std qlen slen' -max_target_seqs 10000 -perc_identity 90 -o <my_blast.tsv> -num_threads 32

#Next, calculate pairwise ANI by combining local alignments between sequence pairs:
python $PHAMB_DIR/checkv/scripts/anicalc.py -i viral.blast.tsv -o ${s}.ani.tsv

#Finally, perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
python $PHAMB_DIR/checkv/scripts/aniclust.py --fna $s.viral.fna --ani ${s}.ani.tsv --out ${s}.viral_cluster.tsv --min_ani 95 --min_tcov 85 --min_qcov 0

#### keep good quality

cut -f 1 $s.viral_cluster.tsv | awk -v s=$s '{printf $0"\t"s"_viral_"; printf("%05d\n", ++i)}' > $workfolder/samples/${s}/virus_detection/$s.viral_renamed.tsv

head -n 1 $workfolder/samples/${s}/virus_detection/unbinned_contigs_new/$s.checkv.out > $s.checkv_oldname.out
cut -f 1 ${s}.viral_cluster.tsv | xargs -I{} grep -w {}  $s.checkv_cand_on.out >> $s.checkv_oldname.out

cut -f 1 ${s}.viral_cluster.tsv | /work_ifs/sukmb276/software/bin/faSomeRecords ${s}.cnd.fna /dev/stdin $s.viral_final.fna

head -n 1 $workfolder/samples/${s}/virus_detection/unbinned_contigs_new/$s.checkv.out >  $s.checkv_final.out
while read line; do
read old new <<< $(echo $line)
grep -w $old $s.checkv_oldname.out | sed "s/$old/$new/" >> $s.checkv_final.out
sed -i -r "s/[>]$old$/\>$new/"  $s.viral_final.fna
done < $workfolder/samples/${s}/virus_detection/$s.viral_renamed.tsv

cp $s.checkv_final.out $s.viral_final.fna ${s}.viral_cluster.tsv $workfolder/samples/${s}/virus_detection/
