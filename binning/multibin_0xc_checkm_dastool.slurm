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

source activate dastool_env
module load ruby

export PATH=$PATH:/work_ifs/sukmb276/software/DAS_Tool/

cd /work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/samples/${s}
mkdir -p dastool
awk '{print $2"\t"$1}' ${s}.maxbin2.contigs_to_bin.tsv > dastool/maxbin2.out
awk '{print $2"\t"$1}' ${s}.metabat2.contigs_to_bin.tsv > dastool/metabat2.out
awk '{print $2"\t"$1}' ${s}.concoct.contigs_to_bin.tsv > dastool/concoct.out
awk '{print $2"\t"$1}' ${s}.vamb.contigs_to_bin.tsv > dastool/vamb.out

cd dastool
cp ../{s}.prodigal.faa .
DAS_Tool -i maxbin2.out,metabat2.out,concoct.out,vamb.out -c $datafolder/Assembly_${s}/${s}.contigs.fa -o ${s}.dastool -p ${s}.prodigal.faa --write_bin_evals --labels=maxbin,metabat,concoct,vamb -t ${SLURM_CPUS_PER_TASK} --write_bins

conda deactivate
source activate binning_env
mkdir -p dastool_bins
cat ${s}.dastool_DASTool_contig2bin.tsv | awk -vs=$s 'BEGIN{i=0; bid=""}{if($2!=bid){i=i+1;bid=$2;nbid=s".dastoolbin_"i".fasta"}; print nbid"\t"$1"\tDASTool"}' > ../${s}.dastool.contig_to_bin.out
cat ../${s}.dastool.contig_to_bin.out | awk '{if(NR==1){print "contig_id,cluster_id"; next}; print $2","$1}' | sed 's/[.]fasta//' | extract_fasta_bins.py $datafolder/Assembly_${s}/${s}.contigs.fa /dev/stdin --output_path dastool_bins

conda deactivate
source load metagenome_env

checkm lineage_wf -x fa --pplacer_threads ${SLURM_CPUS_PER_TASK} -t ${SLURM_CPUS_PER_TASK} -f $workfolder/samples/${s}/${s}.dastool.checkm.out --tab_table dastool_bins checkm_out

cd /work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/samples/${s}
rm -r dastool
