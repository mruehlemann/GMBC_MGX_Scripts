#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=240gb
#SBATCH --output=/home/sukmb276/Isilon/Metagenomes/projects/210916_GMBC_MGX/log/%A_%a.out

workfolder="/home/sukmb276/Isilon/Metagenomes/projects/210916_GMBC_MGX/genome_collection"
cd $workfolder

source activate gtdbtk_1.5.0

cd $TMPDIR

batch=$SLURM_ARRAY_TASK_ID

awk -F '/' '{print $0"\t"$NF}' $workfolder/rep_paths.txt | sed 's/[.]fasta$//' | awk -v batch=$batch '{if(NR%10==batch) print}' > batchfile

GTDBTK_DATA_PATH=/work_ifs/ikmb_repository/shared/microbiome/databases/GTDB/release202/

gtdbtk classify_wf --batchfile batchfile --out_dir core_batch_${batch} --cpus 24

mkdir -p $workfolder/gtdb_tax/core_batch_${batch}

cp core_batch_${batch}/*tsv $workfolder/gtdb_tax/core_batch_${batch}
cp core_batch_${batch}/*user_msa.fasta $workfolder/gtdb_tax/core_batch_${batch}

exit
