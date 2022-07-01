#!/bin/bash
#SBATCH -c 8
#SBATCH --mem=80gb
#SBATCH --time=1-00:00
#SBATCH --output=/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/log/%A_%a.out
#SBATCH --job-name="mb_04b_vamb"

###########################
####    SETUP     #######
###########################
echo $SLURMD_NODENAME

workfolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/"
#cd $workfolder

datafolder="/work_ifs/ikmb_repository/shared/microbiome/GMbC/Assemblies_megahit"
datafolder_reads="/work_ifs/ikmb_repository/shared/microbiome/rawdata/gmbc_data/NovaSeq_Data/"

#source activate metagenome_env

#ls $workfolder/assemblies > $workfolder/samples_available.txt

all_local=($(cut -f 3 /work_ifs/sukmb276/Metagenomes/projects/210916_GMBC_MGX/metadata_gmbc_bn10_complete.tsv | grep -v localit | sort | uniq))

####################################
####    GROUP SELECTION     ########
####################################

mkdir -p $workfolder/allgroups

###########################
####    VAMB        #######
###########################

module load miniconda3
source activate binning_env

head -n 1 $workfolder/samples/$(ls $workfolder/samples | head -n 1)/*.refined.out > $workfolder/allgroups/GMbC.rep97.refined.out
head -n 1 $workfolder/samples/$(ls $workfolder/samples | head -n 1)/*.refined.out > $workfolder/allgroups/GMbC.rep97_mq.refined.out

for this_local in ${all_local[@]};do

cd $workfolder/subgroups/$this_local/

tail -n+2 ${this_local}.dRep_cluster97_representatives.tsv | awk '{if($4>=.7) print $1}'  | xargs -I {} grep {} ${this_local}.cleanbins_list >>  $workfolder/allgroups/GMbC.rep97.cleanbins_list
tail -n+2 ${this_local}.dRep_cluster97_representatives.tsv | awk '{if($4>=.7) print $1}'  | xargs -I {} grep {} ${this_local}.refined.out  >>  $workfolder/allgroups/GMbC.rep97.refined.out
tail -n+2 ${this_local}.dRep_cluster97_representatives.tsv | awk '{if($4>=.5 && $4<.7) print $1}'  | xargs -I {} grep {} ${this_local}.cleanbins_list >>  $workfolder/allgroups/GMbC.rep97_mq.cleanbins_list
tail -n+2 ${this_local}.dRep_cluster97_representatives.tsv | awk '{if($4>=.5 && $4<.7) print $1}'  | xargs -I {} grep {} ${this_local}.refined.out  >>  $workfolder/allgroups/GMbC.rep97_mq.refined.out

done


cd $TMPDIR
cut -f 1 $workfolder/allgroups/GMbC.rep97.cleanbins_list | sed '/^$/d' > genomes_list
dRep compare -g genomes_list -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -sa 0.95 -d .
#cp data_tables/Cdb.csv $workfolder/allgroups/GMbC.rep97.dRep_95.csv
python $workfolder/GMBC_MGX_Scripts/binning/multibin_07a_drep_fastANI99.py
#cp data_tables/Cdb_97.csv $workfolder/allgroups/GMbC.rep97.dRep_97.csv
#cp data_tables/Cdb_99.csv $workfolder/allgroups/GMbC.rep97.dRep_99.csv
#conda deactivate
module load R

R --vanialla <<< '
library(tidyverse)
clusters = read.table("data_tables/Cdb.csv", stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_95=paste0(secondary_cluster)) %>% select(genome, cluster_95) %>% mutate(genome=gsub(".fa$","",genome))

scores = read.table("/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/allgroups/GMbC.rep97.refined.out", head=T, stringsAsFactors=F) %>% tibble::rownames_to_column("genome")  %>% mutate(score.gtdb_rel207_ar53 = ifelse(is.na(score.gtdb_rel207_ar53), 0, score.gtdb_rel207_ar53), score.gtdb_rel207_bac120 = ifelse(is.na(score.gtdb_rel207_bac120), 0, score.gtdb_rel207_ar53)) %>%
	mutate(domain=ifelse(max==score.gtdb_rel207_ar53,"Archaea","Bacteria"), genome=gsub(".fasta","",genome)) %>% select(genome, domain, score=max)

genome_paths=read.table("genomes_list", head=F, stringsAsFactors=F)
genome_paths$genome = sapply(genome_paths$V1, function(x) gsub(".fa","",rev(strsplit(x,split="/")[[1]])[1]))

clusters = left_join(clusters, scores) %>% arrange(-score) %>% distinct_at(.vars=c("cluster_95"), .keep_all=T) %>% left_join(genome_paths)

write.table(clusters$V1, "fastani_ref_list", col.names=F, row.names=F, quote=F)
'

cut -f 1 $workfolder/allgroups/GMbC.rep97_mq.cleanbins_list > fastani_query_list
fastANI --refList fastani_ref_list --ql fastani_query_list -t $SLURM_CPUS_PER_TASK --fragLen 10000 --minFraction 0.8 -o fastani_mq.out

R --vanilla <<< '
library(tidyverse); fa=read.table("fastani_mq.out",head=F, stringsAsFactors=F);
fa$V1=sapply(fa$V1, function(x) rev(strsplit(x,split="/")[[1]])[1]);
fa$V2=sapply(fa$V2, function(x) rev(strsplit(x,split="/")[[1]])[1]);
fa_filt = fa %>% filter(fa$V3 >= 95) ;
mq=read.table("fastani_query_list", head=F, stringsAsFactors=F);
mq$genome = sapply(mq$V1, function(x) rev(strsplit(x,split="/")[[1]])[1]) ;
mq_nohit = mq %>% filter(!genome %in% fa_filt$V1);
write.table(mq_nohit$V1, "drep_l2_input", col.names=F, row.names=F, quote=F)
'

dRep compare -g drep_l2_input -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI -sa 0.95 -d drep_l2_out

R --vanilla <<< '
library(tidyverse);
clusters = read.table("data_tables/Cdb.csv", stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_95=paste0(secondary_cluster)) %>% mutate(genome=gsub(".fa","",genome)) %>% select(genome, cluster_95)

scores = read.table("/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/allgroups/GMbC.rep97.refined.out", head=T, stringsAsFactors=F) %>% tibble::rownames_to_column("genome") %>% mutate(score.gtdb_rel207_ar53 = ifelse(is.na(score.gtdb_rel207_ar53), 0, score.gtdb_rel207_ar53), score.gtdb_rel207_bac120 = ifelse(is.na(score.gtdb_rel207_bac120), 0, score.gtdb_rel207_ar53)) %>%
	mutate(domain=ifelse(max==score.gtdb_rel207_ar53,"Archaea","Bacteria"), genome=gsub(".fasta","",genome)) %>% select(genome, domain, score=max)

clusters = left_join(clusters, scores)

scores2 = read.table("/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/allgroups/GMbC.rep97_mq.refined.out", head=T, stringsAsFactors=F) %>% tibble::rownames_to_column("genome") %>% mutate(score.gtdb_rel207_ar53 = ifelse(is.na(score.gtdb_rel207_ar53), 0, score.gtdb_rel207_ar53), score.gtdb_rel207_bac120 = ifelse(is.na(score.gtdb_rel207_bac120), 0, score.gtdb_rel207_ar53)) %>%
	mutate(domain=ifelse(max==score.gtdb_rel207_ar53,"Archaea","Bacteria"), genome=gsub(".fasta","",genome)) %>% select(genome, domain, score=max)

fa=read.table("fastani_mq.out",head=F, stringsAsFactors=F);
fa$V1=sapply(fa$V1, function(x) gsub(".fa","",rev(strsplit(x,split="/")[[1]])[1]));
fa$V2=sapply(fa$V2, function(x) gsub(".fa","",rev(strsplit(x,split="/")[[1]])[1]));
fa_filt = fa %>% filter(fa$V3 >= 95) %>% left_join(scores2, by=c("V1"="genome")) %>% left_join(clusters %>% select(-score, -domain), by=c("V2"="genome")) %>% arrange(-V3) %>% distinct_at(.vars="V1", .keep_all=T) %>% select(genome=V1, cluster_95, domain, score)

clusters = bind_rows(clusters, fa_filt) %>% arrange(-score)

clusters2 = read.table("drep_l2_out/data_tables/Cdb.csv", stringsAsFactors=F, head=T, sep=",") %>% mutate(cluster_95=paste0(secondary_cluster)) %>% select(genome, cluster_95) %>% mutate(genome=gsub(".fa$","",genome))

clusters2 = left_join(clusters2, scores2) %>% left_join(clusters2 %>% group_by(cluster_95) %>% summarize(n=n())) %>% filter(n>1) %>% mutate(cluster_95=paste0(cluster_95,"_MQ")) %>% arrange(-score)

clusters = bind_rows(clusters, clusters2 %>% select(-n))

cl95_reps = clusters %>% arrange(-score, cluster_95) %>% distinct_at(.vars="cluster_95", .keep_all=T) %>% arrange(-score) %>% mutate(cluster_95_final=paste0("GMbC","_cluster95_",formatC(seq_along(cluster_95), width = 6, format = "d", flag = "0")))

clusters = clusters %>% left_join(cl95_reps %>% select(cluster_95, cluster_95_final))

write.table(clusters, paste0("GMbC", ".dRep_cluster.tsv"), sep="\t", row.names=F,quote=F)
write.table(cl95_reps, paste0("GMbC", ".dRep_cluster95_representatives.tsv"), sep="\t", row.names=F,quote=F)

all_sub = list.files("/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning/subgroups/")
all_sub = all_sub[all_sub!="alltax"]
all_clusters = sapply(all_sub, function(x) read.table(paste0("/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning//subgroups/",x,"/",x,".dRep_cluster.tsv"), head=T, stringsAsFactors=F), simplify=F) %>% do.call("bind_rows",.)

all_reps97 = sapply(all_sub, function(x) read.table(paste0("/work_ifs/ikmb_repository/shared/microbiome/GMbC_binning//subgroups/",x,"/",x,".dRep_cluster97_representatives.tsv"), head=T, stringsAsFactors=F), simplify=F) %>% do.call("rbind",.)

all_clusters_w_rep = all_clusters %>% left_join(all_reps97 %>% select(cluster_97_final, genome_rep97=genome)) %>% mutate(genome=gsub(".fa","",genome), genome_rep97=gsub(".fa","",genome_rep97)) %>% left_join(clusters %>% select(genome_rep97=genome, GMbC_95_final = cluster_95_final))

write.table(all_clusters_w_rep, paste0("GMbC", ".all_genomes.tsv"), sep="\t", row.names=F,quote=F)
'

cp GMbC.* $workfolder/allgroups

cd /$workfolder/allgroups/


while read line; do
read genome cluster_95 domain score cluster_95_final <<< $(echo $line)
grep -w $genome GMbC.rep97.cleanbins_list | awk -v cl=$cluster_95_final '{print $1"\t"cl}'
done < GMbC.dRep_cluster95_representatives.tsv | tee GMbC.representatives_list
