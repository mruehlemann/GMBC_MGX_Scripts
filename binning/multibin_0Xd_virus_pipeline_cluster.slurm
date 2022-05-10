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
####    MAXBIN DEFAULT #######
#############################
cd $TMPDIR

source activate phamb-1.0.2_env
PHAMB_DIR="/work_ifs/ikmb_repository/shared/microbiome/databases/phamb_ressources"

ava=$(ls -d $datafolder/Assembly_${s})
if [ "$ava" == "" ]; then exit; fi
if [ -e "$workfolder/samples/${s}/virus_detection/unbinned_contigs/${s}.checkv.out" ]; then exit; fi

mkdir phamb_discovery

grep '>' $workfolder/samples/${s}/cleanbins/*fa | cut -d '>' -f 2 | cut -d ' ' -f 1 > phamb_discovery/bact_contigs.txt
cat $datafolder/Assembly_${s}/${s}.contigs.fa | cut -d ' ' -f 1 | awk 'BEGIN {RS = ">" ; ORS = ""} length($2) < 2000 {print $1"\n"}' > phamb_discovery/short_contigs.txt
cat phamb_discovery/bact_contigs.txt phamb_discovery/short_contigs.txt | cut -d ' ' -f 1 > phamb_discovery/bad_contigs.txt

/work_ifs/sukmb276/software/bin/faSomeRecords $datafolder/Assembly_${s}/${s}.contigs.fa -exclude phamb_discovery/bad_contigs.txt /dev/stdout | cut -d ' ' -f 1 > phamb_discovery/contigs.fna

grep '>' phamb_discovery/contigs.fna | tr -d '>' | cut -f 1 -d ' ' > phamb_discovery/good_contigs.txt

cat phamb_discovery/bad_contigs.txt | cut -f 1 | awk '{for(i=1;i<=1000;++i) print $1"_"i}' | /work_ifs/sukmb276/software/bin/faSomeRecords $workfolder/samples/${s}/${s}.prodigal.faa -exclude /dev/stdin phamb_discovery/proteins.faa

cd phamb_discovery

mkdir -p annotations
python3 $PHAMB_DIR/DeepVirFinder/dvf.py -i contigs.fna -o DVF -l 2000
mv DVF/contigs.fna_gt2000bp_dvfpred.txt annotations/all.DVF.predictions.txt
#prodigal -i contigs.fna -d genes.fna -a proteins.faa -p meta -g 11

hmmsearch --cpu ${SLURM_CPUS_PER_TASK} -E 1.0e-05 -o output.txt --tblout annotations/all.hmmMiComplete105.tbl $PHAMB_DIR/databases/Bact105.hmm proteins.faa
hmmsearch --cpu ${SLURM_CPUS_PER_TASK} -E 1.0e-05 -o output.txt --tblout annotations/all.hmmVOG.tbl $PHAMB_DIR/databases/AllVOG.hmm proteins.faa

gzip contigs.fna

cat $workfolder/samples/${s}/${s}.concoct.contigs_to_bin.tsv | cut -f 1-2 | sed 's/[.]fasta//' | grep -Fwf good_contigs.txt > concoct_clusters.tsv
cat $workfolder/samples/${s}/${s}.metabat2.contigs_to_bin.tsv | awk '{print $1"\t"$2}' | sed 's/[.]fasta//' | grep -Fwf good_contigs.txt > metabat_clusters.tsv
cat $workfolder/samples/${s}/${s}.vamb.contigs_to_bin.tsv | awk '{print "VAMB"$1"\t"$2}' | sed 's/[.]fa//' | grep -Fwf good_contigs.txt > vamb_clusters.tsv

ln -s $PHAMB_DIR/phamb/workflows/mag_annotation mag_annotation

python $PHAMB_DIR/phamb/workflows/mag_annotation/scripts/run_RF.py contigs.fna.gz metabat_clusters.tsv annotations resultsdir_metabat
python $PHAMB_DIR/phamb/workflows/mag_annotation/scripts/run_RF.py contigs.fna.gz concoct_clusters.tsv annotations resultsdir_concoct
python $PHAMB_DIR/phamb/workflows/mag_annotation/scripts/run_RF.py contigs.fna.gz vamb_clusters.tsv annotations resultsdir_vamb

#ls resultdir
#resultsdir/vambbins_aggregated_annotation.txt
#resultsdir/vambbins_RF_predictions.txt
#resultsdir/vamb_bins #Concatenated predicted viral bins - writes bins in chunks to files so there might be several!

#rm proteins.faa

cat resultsdir_*/vamb_bins/vamb_bins.*.fna > viral_candidates.fna

checkv contamination viral_candidates.fna checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv completeness viral_candidates.fna checkv_out -t ${SLURM_CPUS_PER_TASK} -d $PHAMB_DIR/databases/checkv-db-v1.0/
checkv complete_genomes viral_candidates.fna checkv_out
checkv quality_summary viral_candidates.fna checkv_out

mkdir -p $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output
cp resultsdir_concoct/vambbins_aggregated_annotation.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/concoct_aggregated_annotation.txt
cp resultsdir_concoct/vambbins_RF_predictions.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/concoct_RF_predictions.txt
cp resultsdir_metabat/vambbins_aggregated_annotation.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/metabat_aggregated_annotation.txt
cp resultsdir_metabat/vambbins_RF_predictions.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/metabat_RF_predictions.txt
cp resultsdir_vamb/vambbins_aggregated_annotation.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/vamb_aggregated_annotation.txt
cp resultsdir_vamb/vambbins_RF_predictions.txt $workfolder/samples/${s}/virus_detection/unbinned_contigs/phamb_output/vamb_RF_predictions.txt

cp checkv_out/quality_summary.tsv $workfolder/samples/${s}/virus_detection/unbinned_contigs/${s}.checkv.out
cat checkv_out/*.fna > $workfolder/samples/${s}/virus_detection/unbinned_contigs/${s}.checkv.candidates.fna
