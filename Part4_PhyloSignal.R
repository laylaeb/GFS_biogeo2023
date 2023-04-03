
library(adespatial)
library(picante)
library(phylobase)
library(phylosignal)

ASVs<-read.csv("20221026NOMIS_full_asv_sedimentup_filtered_final.csv", row.names=1)
tree<-read.tree("NOMIS_16S_rooted_fasttree.tree")


ASVs.rand<-ASVs[sample(45737,600),]

bd.ASVs<-beta.div(t(ASVs), method="log.chord",nperm=1)
SCBD<-bd.ASVs$SCBD
tree.rand<-prune.sample(t(ASVs.rand),tree)
dat<-c()
dat$SCBD<-bd.ASVs$SCBD
dat<-as.data.frame(dat)
dat.rand<-dat[rownames(dat) %in% rownames(ASVs.rand),, drop=FALSE]


p4d <- phylo4d(tree.rand, dat.rand)
ps.rand<-phyloSignal(p4d = p4d, method = "all")
ps.rand
pc <- phyloCorrelogram(p4d, n.points=25, ci.bs = 99)
plot(pc)
pc.res<-pc$res
pc.res<-data.frame(pc.res)
points(pc.res$X1, pc.res$X4)
pc$res

plot(pc,add=T)

res<-c()
x<-1
repeat{
ASVs.rand<-ASVs[sample(45737,600),]
tree.rand<-prune.sample(t(ASVs.rand),tree)
dat.rand<-dat[rownames(dat) %in% rownames(ASVs.rand),, drop=FALSE]
p4d <- phylo4d(tree.rand, dat.rand)
ps.rand<-phyloSignal(p4d = p4d, method = "I")
res<-rbind(res,c(ps.rand$stat, ps.rand$pvalue))
x<-x+1
if (x == 200){
  break
}
}



plot(pc)
plot(pc.res$X1, pc.res$X4, type="l", xlab="phylogenetic distance", ylab="correlation", lwd=2, col=rgb(0,0,0,0.55))
x<-1
repeat{
  ASVs.rand<-nomis[sample(45737,600),]
  tree.rand<-prune.sample(t(ASVs.rand),tree)
  dat.rand<-dat[rownames(dat) %in% rownames(ASVs.rand),, drop=FALSE]
  p4d <- phylo4d(tree.rand, dat.rand)
  pc <- phyloCorrelogram(p4d, n.points=25, ci.bs = 99)
  pc.res<-pc$res
  pc.res<-data.frame(pc.res)
  points(pc.res$X1, pc.res$X4, type="l", lwd=2, col=rgb(0,0,0,0.55))
  x<-x+1
  if (x ==25){
    break
  }
}




res2<-data.frame(res)

res2<-res2[res2$I !="NaN",]
colnames(res2)<-c("Moran.I","p-value")
median(as.numeric(res2$Moran.I)); quantile(as.numeric(res2$Moran.I),0.75)

median(as.numeric(res2$`p-value`))
