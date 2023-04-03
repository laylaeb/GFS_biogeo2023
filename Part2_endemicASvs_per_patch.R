
setwd(choose.dir())
ASVs<-read.csv("20221026NOMIS_full_asv_sedimentup_filtered_final.csv", head=T, row.names = 1)
ASVs.unfiltered<-read.csv("20230103_NOMIS_unfiltered_sedimentup_table.csv", sep="",row.names = 1)
tax<-read.csv("20221026NOMIS_full_tax_sedimentup_filtered_final.csv", head=T, row.names=1)
ASVs.endemic<-read.csv("endemic_asv_table.csv")
samp<-read.csv("samples148.csv")





# occurrence of endemic ASVs within replicate patches --------------------

library(reshape2)
#remove samples with < 2 patches (142/148)
X<-table(unlist(strsplit(colnames(ASVs.unfiltered), '_'))[grepl("GL",unlist(strsplit(colnames(ASVs.unfiltered), '_')))])
X1<-X[X>1]
ASVs2<-ASVs[,colnames(ASVs) %in% names(X1)]
write.csv(X,"patches.sampled.csv")


# write .csv files for each region with occurrences of endemic ASVs across patches
for(i in levels(factor(ASVs.endemic$MR))){
  samp.i<-samp[samp$region==i,]
  ASVs.i<-ASVs2[colnames(ASVs2) %in% samp.i$ID]
  ASVs.i.melt<-melt(as.matrix(ASVs.i))
  ASVs.i.melt<-ASVs.i.melt[ASVs.i.melt$value>0,]
  end.i<-ASVs.endemic[ASVs.endemic$MR==i,]
  ASVs.i.occur<-data.frame(matrix(NA, ncol=1, nrow=nrow(end.i)))
  rownames(ASVs.i.occur)<-end.i$ASV
  for(j in levels(factor(ASVs.i.melt$Var2))){
    ASVs.j.melt<-ASVs.i.melt[ASVs.i.melt$Var2==j,]
    ASVs.j.ASV<-ASVs.unfiltered[rownames(ASVs.unfiltered) %in% ASVs.j.melt$Var1,]
    ASVs.j.ASV<-ASVs.j.ASV[,grepl(paste0(j,"_"), colnames(ASVs.j.ASV))]
    ASVs.j.ASV[ASVs.j.ASV>0]<-1
    ASVs.j.occur<-data.frame(rowSums(ASVs.j.ASV))
    colnames(ASVs.j.occur)<-"occur"
    ASVs.i.occur<-merge(ASVs.i.occur,ASVs.j.occur,by="row.names",all.x=TRUE)
    rownames(ASVs.i.occur)<-ASVs.i.occur$Row.names
    ASVs.i.occur<-ASVs.i.occur[,-1]
  }
  ASVs.i.occur<-ASVs.i.occur[,-1]
  colnames(ASVs.i.occur)<-levels(factor(ASVs.i.melt$Var2))
  ASVs.i.occur[is.na(ASVs.i.occur)]<-0
  write.csv(ASVs.i.occur, paste0(i,"_endemics_patches.csv"))
}
