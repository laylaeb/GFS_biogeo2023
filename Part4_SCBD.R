setwd(choose.dir())
ASVs<-read.csv("20221026NOMIS_full_asv_sedimentup_filtered_final.csv", head=T, row.names = 1)
ASVs.unfiltered<-read.csv("20230103_NOMIS_unfiltered_sedimentup_table.csv", sep="",row.names = 1)
tax<-read.csv("20221026NOMIS_full_tax_sedimentup_filtered_final.csv", head=T, row.names=1)
ASVs.endemic<-read.csv("endemic_asv_table.csv")
samp<-read.csv("samples148.csv")
tree<-read.tree("NOMIS_16S_rooted_fasttree.tree")



# taxa contribution to beta-diversity (as total variance) --------------------------------------------------------------------


bd.ASVs<-beta.div(t(ASVs), method="log.chord",nperm=1)

SCBD<-bd.ASVs$SCBD
summary(rownames(tax)==names(SCBD))
tax$SCBD<-SCBD

write.csv(tax, "genus_SCBD.csv")
# use pivot table in Excel to obtain #ASVs, average SCBD per genus and nr of phyloscore <-2 ASVs
gen<-read.csv("genus_ASVs_SCBD.csv", head=T)

gen0<-gen[gen$nr_phyloscore_ASVs==0,]
gen2<-gen[gen$nr_phyloscore_ASVs>0,]
plot(gen$SCBD, gen$ASVs, type="n", log="xy", xlab="SCBD",ylab="number of ASVs per genus")
points(gen0$SCBD, gen0$ASVs, pch=16, col="lightgrey", cex=1.5)
points(gen2$SCBD, gen2$ASVs, pch=21, bg=rgb(50/255,52/255,217/255,gen2$rel_physcore_ASVs), cex=1.5)
text(gen2$SCBD, gen2$ASVs, gen2$genus,cex=0.2 )


