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

all_local=($(cut -f 3 metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

echo $SLURMD_NODENAME
####################################
####    GROUP SELECTION     ########
####################################

this_local=${all_local[$SLURM_ARRAY_TASK_ID]}
mkdir -p $workfolder/subgroups/${this_local}
local_samples=$(awk -veth=${this_local} '{if($3==eth) print $1}' $workfolder/metadata_gmbc_bn10_complete.tsv)
cd $TMPDIR


###########################
####    MAPPING     #######
###########################

mkdir -p bams

for s in ${local_samples}; do
ava=$(grep $s $workfolder/samples_available.txt)
if [ "$ava" == "" ]; then continue; fi
echo $s
minimap2 -t ${SLURM_CPUS_PER_TASK} -N 50 -ax sr  $workfolder/subgroups/${this_local}/${this_local}.catalogue.mmi $datafolder/${s}_P_1.fastq.gz $datafolder/${s}_P_2.fastq.gz > $s.minimap.sam
samtools view -F 3584 -b --threads 20 $s.minimap.sam | samtools sort > bams/$s.minimap.bam
rm $s.minimap.sam
done

jgi_summarize_bam_contig_depths --outputDepth $workfolder/subgroups/${this_local}/${this_local}.allsamples_depth.txt bams/*

### maxbin2 binning GTDB bacteria
mkdir -p $workfolder/subgroups/${this_local}/abufiles
nsmp=$(head -n 1 $workfolder/subgroups/${this_local}/${this_local}.allsamples_depth.txt | awk '{print (NF-3)/2}')
for i in $(seq 1 $nsmp); do
tail -n+2 $workfolder/subgroups/${this_local}/${this_local}.allsamples_depth.txt | awk -vi=$i '{print $1"\t"$(4+(i-1))}' >  $workfolder/subgroups/${this_local}/abufiles/${this_local}.$i.multi.maxbin.cov
done
ls $workfolder/subgroups/${this_local}/abufiles/*.cov > $workfolder/subgroups/${this_local}/${this_local}.abulist
