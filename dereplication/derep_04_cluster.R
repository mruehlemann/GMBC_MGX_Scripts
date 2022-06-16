
.libPaths("/home/sukmb276/.conda/envs/metagenome_env/lib/R/library")

library(ape)
#library(vegan)
library(phangorn)
#library(gtools)

c<-as.matrix(readDist("final.phylip"))
c.sub=c
saveRDS(c, "distmat.Rds")

### initial clustering
hc<-hclust(as.dist(c), "average")

### tree cutting
### values are taken as estimates from GTDB r95
cluster=data.frame(row.names=rownames(c),stl=cutree(hc, h=.01),sl=cutree(hc, h=.065), gl=cutree(hc, h=.23), fl=cutree(hc, h=.375))

cluster_final=cluster
cluster_final[]<-NA
colnames(cluster_final)<-c("StGB","SGB","GGB","FGB")

### load genome stats
checkm<-read.table("mq.checkm.out",head=F, stringsAsFactors=F, sep="\t", row.names=1)
#rownames(checkm)<-gsub("-",".",rownames(checkm))
checkm.sub<-checkm[rownames(c),]

### start with SGB refinement
cl_rep=data.frame(cl=unique(cluster$sl), rep=NA, comp=NA, cont=NA, het=NA)


### selection of best representeative based on checkm statistics
for(u in unique(cluster$sl)){
checkm.sub.sub<-checkm.sub[rownames(cluster[cluster$sl==u,]),]
rep=rownames(checkm.sub.sub)[which.max(rank(checkm.sub.sub$V12) + rank(-checkm.sub.sub$V13) + rank(-checkm.sub.sub$V14))]
cl_rep[cl_rep$cl==u,"rep"]<-rep
cl_rep[cl_rep$cl==u,c("comp","cont","het")]<-checkm.sub.sub[rep,c("V12","V13","V14")]
}


c.sub.rep<-c.sub[cl_rep$rep, cl_rep$rep]
cl_rep$merge=sapply(cl_rep$cl, function(x){cl_new=which(c.sub.rep[x,]<0.05);
        if(length(cl_new)==1){return(cl_new)}; cl_new2=which(apply(c.sub.rep[cl_new,]<0.05,2,any));
        while(!identical(cl_new,cl_new2)){cl_new=cl_new2;cl_new2=which(apply(c.sub.rep[cl_new,]<0.05,2,any))}; cl_new=cl_new2;
        return(paste0(cl_new2,collapse=","))})

cl_rep$SGB=paste0("SGB_",as.numeric(as.factor(cl_rep$merge)))

cluster_final$SGB<-cl_rep[match(cluster$sl, cl_rep$cl),"SGB"]

cl_rep_final=data.frame(cl=unique(cluster_final$SGB), rep=NA, comp=NA, cont=NA, het=NA)

for(u in unique(cluster_final$SGB)){
checkm.sub.sub<-checkm.sub[rownames(cluster_final[cluster_final$SGB==u,]),]
rep=rownames(checkm.sub.sub)[which.max(rank(checkm.sub.sub$V12) + rank(-checkm.sub.sub$V13) + rank(-checkm.sub.sub$V14))]
cl_rep_final[cl_rep_final$cl==u,"rep"]<-rep
cl_rep_final[cl_rep_final$cl==u,c("comp","cont","het")]<-checkm.sub.sub[rep,c("V12","V13","V14")]
}

### GGB
### start with SGB refinement
cl_rep_gl=data.frame(cl=unique(cluster$gl), rep=NA, comp=NA, cont=NA, het=NA)

for(u in unique(cluster$gl)){
checkm.sub.sub<-checkm.sub[rownames(cluster[cluster$gl==u,]),]
rep=rownames(checkm.sub.sub)[which.max(rank(checkm.sub.sub$V12) + rank(-checkm.sub.sub$V13) + rank(-checkm.sub.sub$V14))]
cl_rep_gl[cl_rep_gl$cl==u,"rep"]<-rep
cl_rep_gl[cl_rep_gl$cl==u,c("comp","cont","het")]<-checkm.sub.sub[rep,c("V12","V13","V14")]
}


#c.sub.rep_gl<-c.sub[cl_rep_gl$rep, cl_rep_gl$rep]

#cl_rep_gl$merge=sapply(cl_rep_gl$cl, function(x){cl_new=which(c.sub.rep_gl[x,]<0.23);
#        if(length(cl_new)==1){return(cl_new)}; cl_new2=which(apply(c.sub.rep_gl[cl_new,]<0.23,2,any));
#        while(!identical(cl_new,cl_new2)){cl_new=cl_new2;cl_new2=which(apply(c.sub.rep_gl[cl_new,]<0.23,2,any))}; cl_new=cl_new2;
#        return(paste0(cl_new2,collapse=","))})

cl_rep_gl$GGB=paste0("GGB_",as.numeric(as.factor(cl_rep_gl$cl)))

cluster_final$GGB<-cl_rep_gl[match(cluster$gl, cl_rep_gl$cl),"GGB"]

cl_rep_final_gl=data.frame(cl=unique(cluster_final$GGB), rep=NA, comp=NA, cont=NA, het=NA)

for(u in unique(cluster_final$GGB)){
checkm.sub.sub<-checkm.sub[rownames(cluster_final[cluster_final$GGB==u,]),]
rep=rownames(checkm.sub.sub)[which.max(rank(checkm.sub.sub$V12) + rank(-checkm.sub.sub$V13) + rank(-checkm.sub.sub$V14))]
cl_rep_final_gl[cl_rep_final_gl$cl==u,"rep"]<-rep
cl_rep_final_gl[cl_rep_final_gl$cl==u,c("comp","cont","het")]<-checkm.sub.sub[rep,c("V12","V13","V14")]
}


### FGB
### start with SGB refinement
cl_rep_fl=data.frame(cl=unique(cluster$fl), rep=NA, comp=NA, cont=NA, het=NA)

for(u in unique(cluster$fl)){
checkm.sub.sub<-checkm.sub[rownames(cluster[cluster$fl==u,]),]
rep=rownames(checkm.sub.sub)[which.max(rank(checkm.sub.sub$V12) + rank(-checkm.sub.sub$V13) + rank(-checkm.sub.sub$V14))]
cl_rep_fl[cl_rep_fl$cl==u,"rep"]<-rep
cl_rep_fl[cl_rep_fl$cl==u,c("comp","cont","het")]<-checkm.sub.sub[rep,c("V12","V13","V14")]
}


#c.sub.rep_fl<-c.sub[cl_rep_fl$rep, cl_rep_fl$rep]
#cl_rep_fl$merge=sapply(cl_rep_fl$cl, function(x){cl_new=which(c.sub.rep_fl[x,]<0.3);
#        if(length(cl_new)==1){return(cl_new)}; cl_new2=which(apply(c.sub.rep_fl[cl_new,]<0.3,2,any));
#        while(!identical(cl_new,cl_new2)){cl_new=cl_new2;cl_new2=which(apply(c.sub.rep_fl[cl_new,]<0.3,2,any))}; cl_new=cl_new2;
#        return(paste0(cl_new2,collapse=","))})

cl_rep_fl$FGB=paste0("FGB_",as.numeric(as.factor(cl_rep_fl$cl)))

cluster_final$FGB<-cl_rep_fl[match(cluster$fl, cl_rep_fl$cl),"FGB"]

cl_rep_final_fl=data.frame(cl=unique(cluster_final$FGB), rep=NA, comp=NA, cont=NA, het=NA)

for(u in unique(cluster_final$FGB)){
checkm.sub.sub<-checkm.sub[rownames(cluster_final[cluster_final$FGB==u,]),]
rep=rownames(checkm.sub.sub)[which.max(rank(checkm.sub.sub$V12) + rank(-checkm.sub.sub$V13) + rank(-checkm.sub.sub$V14))]
cl_rep_final_fl[cl_rep_final_fl$cl==u,"rep"]<-rep
cl_rep_final_fl[cl_rep_final_fl$cl==u,c("comp","cont","het")]<-checkm.sub.sub[rep,c("V12","V13","V14")]
}
####

write.table(cl_rep_final_fl,"FGB_representatives.tsv", sep="\t", quote=F)
write.table(cl_rep_final_gl,"GGB_representatives.tsv", sep="\t", quote=F)
write.table(cl_rep_final,"SGB_representatives.tsv", sep="\t", quote=F)


write.table(cluster_final, "fastANI_cluster_AvLinkage_refined.tsv", sep="\t", quote=F)
