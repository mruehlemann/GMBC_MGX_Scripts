#!/bin/bash
#SBATCH -c 4
#SBATCH --mem=40gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_03a_cds"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

source activate metagenome_env

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
cd $TMPDIR

######################################################
######### CDS CALLING AND ANNOTATION     ############
######################################################
#for s in ${local_samples}; do
#ava=$(grep $s $workfolder/samples_available.txt)
ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/${s}.hmm" ]; then exit; fi
echo $s
mkdir -p $workfolder/samples/${s}
mkdir -p ${s}_tmp_workfolder
cat $datafolder/Assembly_${s}/${s}.contigs.fa | ~/Isilon/software/bin/parallel -j ${SLURM_CPUS_PER_TASK} --block 999k --recstart '>' --pipe prodigal -p meta -a ${s}_tmp_workfolder/${s}.{#}.faa -d ${s}_tmp_workfolder/${s}.{#}.ffn -o tmpfile
cat ${s}_tmp_workfolder/$s.*.faa > $workfolder/samples/${s}/${s}.prodigal.faa
cat ${s}_tmp_workfolder/$s.*.ffn > $workfolder/samples/${s}/${s}.prodigal.ffn

hmmsearch -o ${s}.hmm.tigr.out --tblout ${s}.hmm.tigr.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/tigrfam/tigrfam.hmm $workfolder/samples/${s}/${s}.prodigal.faa
#pfam_scan.pl -cpu ${SLURM_CPUS_PER_TASK} -outfile ${sample}.hmm.pfam.hit.out -fasta $outdir/$sample/${sample}.prodigal.faa -dir /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/
hmmsearch -o ${s}.hmm.pfam.out --tblout ${s}.hmm.pfam.hit.out --noali --notextw --cut_nc --cpu ${SLURM_CPUS_PER_TASK} /work_ifs/ikmb_repository/databases/GTDB-TK/release207/markers/pfam/Pfam-A.hmm $workfolder/samples/${s}/${s}.prodigal.faa

cat ${s}.hmm.tigr.hit.out | grep -v "^#" | awk '{print $1"\t"$3"\t"$5}' > ${s}.tigr
cat ${s}.hmm.pfam.hit.out | grep -v "^#" | awk '{print $1"\t"$4"\t"$5}' > ${s}.pfam
cat ${s}.pfam ${s}.tigr > $workfolder/samples/${s}/${s}.hmm

jgi_summarize_bam_contig_depths --outputDepth $workfolder/samples/${s}/${s}.depth.txt $datafolder/Assembly_${s}/$s.sorted.bam
