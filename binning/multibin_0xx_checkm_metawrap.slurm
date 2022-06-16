#!/bin/bash
#SBATCH -c 6
#SBATCH --mem=60Gb
#SBATCH --time=10-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_06x_checkm"

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

this_local=$(grep -w $s /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | cut -f 3 )

#############################
####    MAXBIN DEFAULT #######
#############################
cd $TMPDIR

ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/dastool_out/checkm.out" ]; then exit; fi

cat $datafolder/Assembly_${s}/${s}.contigs.fa | cut -d ' ' -f 1 > contigs.fa

source activate binning_env
mkdir -p metabat2_bins
cat $workfolder/samples/${s}/${s}.metabat2.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py contigs.fa /dev/stdin --output_path metabat2_bins
mkdir -p maxbin2_bins
cat $workfolder/samples/${s}/${s}.maxbin2.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py contigs.fa /dev/stdin --output_path maxbin2_bins
mkdir -p concoct_bins
cat $workfolder/samples/${s}/${s}.concoct.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py contigs.fa /dev/stdin --output_path concoct_bins
mkdir -p vamb_bins
cat $workfolder/samples/${s}/${s}.vamb.contigs_to_bin.tsv | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py contigs.fa /dev/stdin --output_path vamb_bins

source activate metawrap_env
export PATH=$PATH:/work_ifs/sukmb276/github/metaWRAP/bin

/work_ifs/sukmb276/github/metaWRAP/bin/metawrap-modules/bin_refinement2.sh -o vamb_and_maxbin -A vamb_bins -B maxbin2_bins -t 24 -m 240 -c 50
/work_ifs/sukmb276/github/metaWRAP/bin/metawrap-modules/bin_refinement2.sh -o metabat_and_concoct -A metabat2_bins -B concoct_bins -t 24 -m 240 -c 50
/work_ifs/sukmb276/github/metaWRAP/bin/metawrap-modules/bin_refinement2.sh -o all -A vamb_and_maxbin/metawrap_50_10_bins -B metabat_and_concoct/metawrap_50_10_bins -t 24 -m 240 -c 50

cat all/metawrap_50_10_bins.contigs | awk -v s=$s '{print "MetaWRAP_"s"_"$2"\t"$1"\tmetawrap"}' >  $workfolder/samples/${s}/${s}.metawrap.contigs_to_bin.tsv
