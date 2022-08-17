library(tidyverse)

allsamples = list.files("samples")

alltax = lapply(allsamples, function(x) read.table(paste0("samples/",x,"/",x,".final_tax.tsv"), head=T, stringsAsFactors=F)) %>% do.call("bind_rows",.) %>% data.frame(stringsAsFactors=F)

dir.create("allgroups/tax_tables/", recursive=T)

for(lv in unique(alltax$level)){
	print(lv)
	alltax %>% filter(level == lv) %>% reshape2::dcast(sample ~ tax, value.var="tpm", fill=0) %>% data.frame(row.names=1) %>% write.table(., paste0("allgroups/tax_tables/GMbC_final_table_",lv,".tsv"), sep="\t")
}


allsmry = lapply(allsamples, function(x) read.table(paste0("samples/",x,"/",x,".genome_summary.tsv"), head=T, stringsAsFactors=F)) %>% do.call("bind_rows",.) %>% data.frame(stringsAsFactors=F)
write.table(allsmry, "allgroups/GMbC.all_genome_summary.tsv", sep="\t",row.names=F)
