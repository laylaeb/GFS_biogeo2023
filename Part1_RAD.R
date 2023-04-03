library(vegan)
library(stringr)

ASVs<-read.csv("20221026NOMIS_full_asv_sedimentup_filtered_final.csv", row.names=1)
tax<-read.csv("20221026NOMIS_full_tax_sedimentup_filtered_final.csv", row.names=1)

ASVs.rel<-decostand(ASVs, "total",MARGIN=2)
colSums(ASVs.rel)

# filtered dataset
par(mfrow=c(12,13), mar=c(0,0.5,1,0),oma=c(0.1,0.1,0.1,0.1), cex.main=0.75)

ASVs.rel<-ASVs.rel[,str_sort(colnames(ASVs.rel), numeric = TRUE)]
for(i in 1:ncol(ASVs.rel)){
  ASVs.i<-ASVs.rel[,i]
  rad.ln.i<-rad.lognormal(ASVs.i, family=gaussian)
  rad.null.i<-rad.null(ASVs.i)
  plot(rad.ln.i, col="grey", main=colnames(ASVs.rel)[i], xaxt="n", yaxt="n") 
  lines(rad.ln.i, col="red", lwd=2)
  lines(rad.null.i, col="cyan", lwd=2)
}



#unfiltered dataset
ASVs.unfilt<-read.csv("20230103_NOMIS_unfiltered_sedimentup_table.csv",sep="", head=T,row.names = 1)
ASVs.unfilt<-ASVs.unfilt[,str_sort(colnames(ASVs.unfilt), numeric = TRUE)]

colnames(ASVs.unfilt)<-sub("\\_.*", "", colnames(ASVs.unfilt))


nms<-unique(colnames(ASVs.unfilt))

par(mfrow=c(12,13), mar=c(0,0.5,1,0),oma=c(0.1,0.1,0.1,0.1), cex.main=0.75)
for(i in 1:148){
  nms.i<-nms[i]
  ASVs.i<-data.frame(ASVs.unfilt[,colnames(ASVs.unfilt) %in% nms.i])
  ASVs.i<-rowMeans(ASVs.i)
  rad.ln<-rad.lognormal(ASVs.i, family=gaussian)
  plot(rad.ln, col="grey", main=nms.i, xaxt="n", yaxt="n") 
  rad.null<-rad.null(ASVs.i)
  lines(rad.ln, col="red", lwd=2)
  lines(rad.null, col="cyan", lwd=2)
}

