library(tidyverse)
reps=read.table("SGB_representatives.tsv", head=T, stringsAsFactors=F)
reps$cl = paste0("GMBC_", reps$cl)

tax=read.table("gtdb_tax/alltax.summary.tsv", head=T, stringsAsFactors=F, sep="\t")

pas<-read.table("pasolli_mash.out", head=F, stringsAsFactors=F)
pas$V2 = gsub("[.]fasta","", pas$V2)
pas$V2 = gsub("mash_comp_ref/","", pas$V2)
pas$V1 = gsub("representatives/|[.]fa","", pas$V1)
pas_sub = pas %>% arrange(V3) %>% filter(!duplicated(V2)) %>% dplyr::select(V1,V2,V3)

man<-read.table("manara_mash.out", head=F, stringsAsFactors=F)
man$V2 = gsub("mash_comp_ref/","", man$V2)
man$V2 = gsub("[.]fasta","", man$V2)
man$V1 = gsub("[.]fa","", man$V1)
man_sub = man %>% arrange(V3) %>% filter(!duplicated(V2)) %>% dplyr::select(V1,V2,V3)

gra<-read.table("greatapes_mash.out", head=F, stringsAsFactors=F)
gra$V2 = gsub("mash_comp_ref/","", gra$V2)
gra$V2 = gsub("[.]fasta","", gra$V2)
gra$V1 = paste0("GreatApe_",gsub("MASH_REF/|[.]fasta","", gra$V1))
gra_sub = gra %>% arrange(V3) %>% filter(!duplicated(V2)) %>% dplyr::select(V1,V2,V3)

reps2 = left_join(reps, tax %>% dplyr::select(user_genome, classification), by=c("rep"="user_genome")) %>% left_join(., gra_sub, by=c("rep"="V2")) %>% left_join(., man_sub, by=c("rep"="V2")) %>% left_join(., pas_sub, by=c("rep"="V2"))
colnames(reps2)[7:12]=c("nearest_GreatApe", "dist_GreatApe", "nearest_Manara", "dist_Manara", "nearest_Passolli", "dist_Passolli")

write.table(reps2, "GMBC_SGB_annot.tsv", quote=F, sep="\t")

cluster=read.table("fastANI_cluster_AvLinkage_refined.tsv", head=T, stringsAsFactors=F)

reps3= rownames_to_column(cluster, "bin") %>% mutate(locality=gsub("_cleanbin_[0-9]+", "", bin)) %>% mutate(SGB=paste0("GMBC_",SGB)) %>% reshape2::dcast(SGB ~ locality, data=.) %>% left_join(reps2, ., by=c("cl" = "SGB"))

write.table(reps3, "GMBC_SGB_annot_locality.tsv", quote=F, sep="\t")
