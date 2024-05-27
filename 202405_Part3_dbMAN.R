ASVs<-read.csv("20230825_NOMIS_rarefied_deblur_table.csv", head=T, row.names = 1)
samp<-read.csv("samples.csv")

# indicator ASVs ----------------------------------------------------------
library(foreach)
library(doParallel)
library(adiv)


#remove Uganda
samp2<-samp[samp$region !="Uganda",]
asv_table<-ASVs[colnames(ASVs) %in% samp2$ID]
asv_table<-asv_table[rowSums(asv_table)>0,]
asv_table<-asv_table[rowSums(asv_table>0)>1,] #remove unique ASVs

n.cores<-detectCores()
registerDoParallel(cl <- makeCluster(n.cores-1))

res.f.all<-foreach(f = levels(factor(metadata$region)), .packages = c("adiv")) %dopar% {
  metadata.f<-metadata
  metadata.f$region[!grepl(f, metadata.f$region)]<-"other"
  Q.f<-dbMANOVAspecies(t(asv_table),metadata.f$region, nrep=999, global=FALSE,species=FALSE, padj="BH",method="BrayCurtis")
  Q.f.pairwise_adj <- dbMANOVAspecies_pairwise(Q.f)
  list(summary = summary(Q.f.pairwise_adj), pairwise = Q.f.pairwise_adj)
}
save.image("db.MANOVA.res.RData")



