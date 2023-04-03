setwd(choose.dir())
ASVs<-read.csv("20221026NOMIS_full_asv_sedimentup_filtered_final.csv", head=T, row.names = 1)
ASVs.unfiltered<-read.csv("20230103_NOMIS_unfiltered_sedimentup_table.csv", sep="",row.names = 1)
tax<-read.csv("20221026NOMIS_full_tax_sedimentup_filtered_final.csv", head=T, row.names=1)
ASVs.endemic<-read.csv("endemic_asv_table.csv")
samp<-read.csv("samples148.csv")
tree<-read.tree("NOMIS_16S_rooted_fasttree.tree")



# phylogenetic structure of beta-diversity --------------------------------

library(picante)
library(adespatial)
library(speedyseq)
library(phyloseq)




# mean ASV abundance per region 
regs<-c()
for(j in levels(factor(samp$region))){
  reg.j<-samp[samp$region==j,]
  n.j<-rowMeans(ASVs[reg.j$ID])
  regs<-cbind(regs,n.j)
}
regs<-data.frame(regs)
colnames(regs)<-levels(factor(samp$region))

regs<-regs[,1:10] #drop Uganda


#Mean nearest taxon distance
MNTDs<-c()
for(k in levels(factor(samp$region))){
  reg.k<-samp[samp$region==k,]
  n.k<-ASVs[reg.k$ID]
  n.k<-n.k[rowSums(n.k)>0,]
  tree.k<-prune.sample(t(n.k), tree)
  mntd.k<-mntd(t(n.k), cophenetic(tree.k))
  mntd.k<-data.frame(mntd.k)
  mntd.k$reg<-k
  MNTDs<-rbind(MNTDs, mntd.k)
}
boxplot(MNTDs$mntd.k~MNTDs$reg, las=2, xlab="", ylab="MNTD")
abline(h=mean(na.omit(as.numeric(MNTDs$mntd.k))), col="red", lwd=2)


# beta diversity profiles across regions 
across.part<-c()
b.part.across.all<-beta.div.comp(t(regs), coef="S", quant=T)
across.part<-rbind(across.part, c(b.part.across.all$part, nrow(regs),0))

regs.phy<-phyloseq(otu_table(regs, taxa_are_rows = TRUE), phy_tree(tree))
for(h in seq(0.005,0.2,0.005)){
  glom.h<-speedyseq::tree_glom(regs.phy, h)
  OTU.h<-as.data.frame(otu_table(glom.h), "matrix")
  b.part.across<-beta.div.comp(t(OTU.h), coef="S", quant=T)
  across.part<-rbind(across.part, c(b.part.across$part, nrow(OTU.h), h))
}
across.part<-data.frame(across.part)
colnames(across.part)[6:7]<-c("tips","phy.agglom")

plot(across.part$phy.agglom, across.part$BDtotal,type="b",lwd=1.75, cex=1, xlim=c(0,0.2), ylim=c(0,0.5), xlab="phylogenetic agglomeration depth", ylab="beta diversity")
abline(v=mean(na.omit(MNTDs$mntd.k)), lty=1, lwd=2)


# beta diversity profiles across regions 

samp2<-samp[samp$region!="Uganda",]
regs.parts<-c()
for(a in levels(factor(samp2$region))){
  samp.a<-samp2[samp2$region==a,]
  reg.a<-ASVs[samp.a$ID]
  reg.a<-reg.a[rowSums(reg.a)>0,]
  b.part.a<-beta.div.comp(t(reg.a), coef="S", quant=T)
  regs.parts<-rbind(regs.parts, c(b.part.a$part, nrow(reg.a), nrow(reg.a)/nrow(reg.a), h=0, a))}
for(a in levels(factor(samp2$region))){
  samp.a<-samp2[samp2$region==a,]
  reg.a<-ASVs[samp.a$ID]
  reg.a<-reg.a[rowSums(reg.a)>0,]
  reg.a.tree<-prune.sample(t(reg.a),tree)
  reg.a.phy<-phyloseq(otu_table(reg.a, taxa_are_rows = TRUE), phy_tree(reg.a.tree))
  for(h in seq(0.005,0.2,0.005)){
    glom.h<-speedyseq::tree_glom(reg.a.phy, h)
    OTU.h<-as.data.frame(otu_table(glom.h), "matrix")
    b.part<-beta.div.comp(t(OTU.h), coef="S", quant=T)
    regs.parts<-rbind(regs.parts, c(b.part$part, nrow(OTU.h), nrow(OTU.h)/nrow(reg.a), h, a))
    print(c(a,h))}
}

regs.parts<-data.frame(regs.parts)
colnames(regs.parts)[6:9]<-c("n_tips","perc_tips","phy_glom","region")


plot(regs.parts$phy_glom,regs.parts$BDtotal, type="n", ylim=c(0,6), xlab="phylogenetic agglomeration depth", ylab="beta diversity change (%)")
for(i in levels(factor(regs.parts$region))){
  reg.i<-regs.parts[regs.parts$region==i,]
  points(reg.i$phy_glom[1:40], abs(diff(as.numeric(reg.i$BDtotal)))/as.numeric(reg.i$BDtotal)[1]*100, pch=16, cex=1, col="grey")
}
points(across.part$phy.agglom[1:40], abs(diff(across.part$BDtotal)/across.part$BDtotal[1])*100, pch=16, col="red", type="b")
abline(v=0.07, lty=2)


x<-as.numeric(across.part$phy.agglom)[1:40]
y<-abs(diff(across.part$BDtotal)/across.part$BDtotal[1])*100

m2<-nls(y~a*exp(b*x),start=list(a=8,b=-50))
lines(x,predict(m2),col="red",lwd=3)

