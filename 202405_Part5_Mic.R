
ASVs<-read.csv("20230825_NOMIS_rarefied_deblur_table.csv", row.names=1)
tree<-read.tree("20230825_NOMIS_rarefied_deblur.tree")
tax<-read.csv("20230825_NOMIS_rarefied_deblur_taxonomy.csv", row.names = 1)
samp<-read.csv("samples.csv")

library(picante)

levels(factor(tax$Genus))
goi<-data.frame(table(tax$Genus))
goi<-goi[order(goi$Freq, decreasing=TRUE), ]

goi2<-goi[3:32,]
goi2<-goi2[!goi2$Var1 %in% c(" g__OM190"," g__vadinHA49"," g__A21b"," g__0319-6G20", " g__TRA3-20"," g__AKYH767", " g__Anaeromyxobacter"," g__Blfdi19"," g__Candidatus_Udaeobacter"),]
goi3<-goi[3:102,]


res<-c()
for(g in levels(factor(goi3$Var1))){
  tax.g<-tax[tax$Genus==g,]
  com.g<-ASVs[rownames(ASVs) %in% rownames(tax.g),]
  com.g<-com.g[,colSums(com.g>0)>2]
  tree.g<-prune.sample(t(com.g),tree)
  coph.g<-mean(cophenetic(tree.g))
  res1<-rbind(res1,c(coph.g,g))}


library(ggplot2)
library(ggridges)

ggplot(res,aes(x=phy.dist, 
                y=fct_reorder(genus,order), 
                fill=genus))+
  geom_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")



#NMDS of microdiverse genera

pol<-tax[tax$Genus==" g__Polaromonas",]
pols<-ASVs[rownames(ASVs) %in% rownames(pol),]
pol.MDS<-metaMDS(t(pols))
pol.anosim<-anosim(t(pols),samp$region)

rhod<-tax[tax$Genus==" g__Rhodoferax",]
rhods<-ASVs[rownames(ASVs) %in% rownames(rhod),]
rhod.MDS<-metaMDS(t(rhods))
rhod.anosim<-anosim(t(rhods),samp$region)

meth<-tax[tax$Genus==" g__Methylotenera",]
meths<-ASVs[rownames(ASVs) %in% rownames(meth),]
meths<-meths[,!(names(meths) %in% "GL146")]
meth.MDS<-metaMDS(t(meths))
meth.samp<-samp[samp$ID!="GL146",]
meth.anosim<-anosim(t(meths),meth.samp$region)

rhizo<-tax[tax$Genus==" g__Rhizobacter",]
rhizo<-ASVs[rownames(ASVs) %in% rownames(rhizo),]
rhizo.MDS<-metaMDS(t(rhizo))
rhizo.anosim<-anosim(t(rhizo),samp$region)

par(mfrow=c(2,2))
plot(pol.MDS, main="Polaromonas")
ordispider(pol.MDS, samp$region)
legend("topright",paste0("Anosim R: ",round(pol.anosim$statistic,2),"\n", "p-value: ",round(pol.anosim$signif,3)), bty="n")
plot(rhod.MDS, main="Rhodoferax")
ordispider(rhod.MDS, samp$region)
legend("topright",paste0("Anosim R: ",round(rhod.anosim$statistic,2),"\n", "p-value: ",round(rhod.anosim$signif,3)), bty="n")
plot(meth.MDS, main="Methylotenera")
ordispider(meth.MDS, meth.samp$region)
legend("topright",paste0("Anosim R: ",round(meth.anosim$statistic,2),"\n", "p-value: ",round(meth.anosim$signif,3)), bty="n")
plot(rhizo.MDS,  main="Rhizobacter")
ordispider(rhizo.MDS, samp$region)
legend("topright",paste0("Anosim R: ",round(rhizo.anosim$statistic,2),"\n", "p-value: ",round(rhizo.anosim$signif,3)), bty="n")



