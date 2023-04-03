setwd(choose.dir())
ASVs<-read.csv("20221026NOMIS_full_asv_sedimentup_filtered_final.csv", head=T, row.names = 1)
ASVs.unfiltered<-read.csv("20230103_NOMIS_unfiltered_sedimentup_table.csv", sep="",row.names = 1)
tax<-read.csv("20221026NOMIS_full_tax_sedimentup_filtered_final.csv", head=T, row.names=1)
ASVs.endemic<-read.csv("endemic_asv_table.csv")
samp<-read.csv("samples148.csv")




# indicator ASVs ----------------------------------------------------------

library(foreach)
library(doParallel)


n.cores<-detectCores()
registerDoParallel(cl <- makeCluster(n.cores-1))


uganda="Uganda"
prune_Uganda <- subset_samples(NOMIS_FR, !Site_c %in% uganda)
table_taxa <- tax_table(prune_Uganda)

asv_table<- as.data.frame((otu_table(prune_Uganda, taxa_are_rows=T)))
asv_table<- asv_table[rowSums(asv_table[])>0,]
metadata<- as.data.frame(sample_data(prune_Uganda))


res.f.all<-foreach(f = levels(factor(metadata$Site_c)), .packages = c("adiv")) %dopar% {
  metadata.f<-metadata
  metadata.f$Site_c[!grepl(f, metadata.f$Site_c)]<-"other"
  Q.f<-dbMANOVAspecies(t(asv_table),metadata.f$Site_c, nrep=999, global=FALSE,species=FALSE, padj="BH",method="BrayCurtis")
  Q.f.pairwise_adj <- dbMANOVAspecies_pairwise(Q.f)
  c(summary(Q.f.pairwise_adj))
}
