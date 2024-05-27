setwd("/Users/hpeter.INTRANET/Desktop/nomis/Leila/nature_revision/phybeta/")

ASVs<-read.csv("20230825_NOMIS_rarefied_deblur_table.csv", row.names=1)
tree<-read.tree("20230825_NOMIS_rarefied_deblur.tree")
tax<-read.csv("20230825_NOMIS_rarefied_deblur_taxonomy.csv", row.names = 1)
samp<-read.csv("samples.csv")

library(picante)

levels(factor(tax$Genus))
goi<-data.frame(table(tax$Genus))
goi<-goi[order(goi$Freq, decreasing=TRUE), ]
head(goi)

goi2<-goi[3:32,]
goi2<-goi2[!goi2$Var1 %in% c(" g__OM190"," g__vadinHA49"," g__A21b"," g__0319-6G20", " g__TRA3-20"," g__AKYH767", " g__Anaeromyxobacter"," g__Blfdi19"," g__Candidatus_Udaeobacter"),]
goi3<-goi[3:102,]


res1<-c()
for(g in levels(factor(goi3$Var1))){
  tax.g<-tax[tax$Genus==g,]
  com.g<-ASVs[rownames(ASVs) %in% rownames(tax.g),]
  com.g<-com.g[,colSums(com.g>0)>2]
  tree.g<-prune.sample(t(com.g),tree)
  coph.g<-mean(cophenetic(tree.g))
  res1<-rbind(res1,c(coph.g,g))}

write.csv(res1,"overall_phy_dists_genera.csv")


res<-c()
for(g in levels(factor(goi2$Var1))){
  tax.g<-tax[tax$Genus==g,]
  com.g<-ASVs[rownames(ASVs) %in% rownames(tax.g),]
  com.g<-com.g[,colSums(com.g>0)>2]
  tree.g<-prune.sample(t(com.g),tree)
  for(k in 1:ncol(com.g)){
    com.k<-com.g[,k, drop=FALSE]
    com.k<-com.k[com.k>0,,drop=FALSE]
    tree.k<-prune.sample(t(com.k),tree.g)
    coph.k<-mean(cophenetic(tree.k))
    res<-rbind(res, c(g, coph.k, colnames(com.k)))}
   }

res<-data.frame(res)
colnames(res)<-c("genus","phy.dist","sample")
par(mar=c(4,12,2,2))
boxplot(as.numeric(res$phy.dist)~res$genus, las=2, horizontal = TRUE, xlab="", ylab="")

write.csv(res, "phylogenetic_distances_genera.csv")

res2<-read.csv("phylogenetic_distances_genera.csv")

library(ggplot2)
library(ggridges)

ggplot(res2,aes(x=phy.dist, 
                y=fct_reorder(genus,order), 
                fill=genus))+
  geom_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")



ggplot(res2, aes(x = as.numeric(phy.dist), y = factor(genus, ordered = T), fill = genus)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")




ggplot(res2,aes(x=phy.dist, 
           y=fct_reorder(genus,order), 
           fill=genus))+
  geom_density_ridges()+
  theme_ridges() + 
  theme(legend.position = "none")


ggplot(res2, aes(x = phy.dist,
         y = fct_reorder(genus, phy.dist, .fun = mean),
         height=,
         fill = fct_reorder(genus, phy.dist, .fun = mean)
       )) + 
  geom_density_ridges(stat="identity",scale = 3) +
  theme(legend.position = "none")






all.anosim<-anosim(t(ASVs), samp$region)


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




res3<-c()
for(g in levels(factor(goi3$Var1))){
  tax.g<-tax[tax$Genus==g,]
  com.g<-ASVs[rownames(ASVs) %in% rownames(tax.g),]
  com.g<-com.g[,colSums(com.g>0)>2]
  samp.g<-samp[samp$ID %in% names(com.g),]
  anosim.g<-anosim(t(com.g),samp.g$region)
  res3<-rbind(res3, c(anosim.g$statistic, anosim.g$signif,g))
    }

res3<-data.frame(res3)
plot(res3$X1)





#NRI vs altitude

library(picante)
library(iCAMP)
library(ggplot2)
library(ggplot2)
library(ggridges)
phydist<-pdist.big(tree, tree.asbig = FALSE, output = FALSE, nworker = 4)


NRI.res<-c()
for(g in levels(factor(goi2$Var1))){
  tax.g<-tax[tax$Genus==g,]
  com.g<-ASVs[rownames(ASVs) %in% rownames(tax.g),]
  com.g<-com.g[,colSums(com.g>0)>2]
  tree.g<-prune.sample(t(com.g),tree)
  coph.g<-cophenetic.phylo(tree.g)
  NRI.k<-ses.mpd(t(com.g), coph.g, 
            null.model = "taxa.labels")
  NRI.res<-rbind(NRI.res, cbind(NRI.k, g ))}


ggplot(NRI.res, aes(x = as.numeric(mpd.obs.z), y = g, fill = g)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

NRI.res.sel<-NRI.res[NRI.res$g %in% c(" g__Polaromonas"," g__Rhizobacter"," g__Methylotenera"," g__Rhodoferax"),]

ggplot(NRI.res.sel, aes(x = as.numeric(mpd.obs.z), y = g, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01) +
  scale_fill_viridis_c(name = "NRI", option = "C")
