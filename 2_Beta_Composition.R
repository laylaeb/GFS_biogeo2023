## Beta diversity - Community Composition - NMDS and statistics
library(vegan)
library(ggplot2)
library(betadisper)
library(adonis)
library(ecodist)
library(phyloseq)

vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

metadata_nmds <- sample.data.frame(prune_Uganda)

metadata_nmds$lat_attribute <- "A"
metadata_nmds$lat_attribute <- ifelse((metadata_nmds$lat_sp) < 1, metadata_nmds$lat_attribute == "B", metadata_nmds$lat_attribute == "A")
metadata_nmds$lat_attribute <- as.factor(metadata_nmds$lat_attribute)

##NMDS plot
asv_table_nmds <- as.matrix((otu_table(prune_Uganda, taxa_are_rows=T)))
asv_table_nmds_f <- asv_table_nmds[rowSums(asv_table_nmds[])>0,]
nmds_bc_nomis <- metaMDS(t(log1p(asv_table_nmds_f)), distance = "bray", k = 2, trymax=999)

stressplot(nmds_bc_nomis)
data.scores = as.data.frame(scores(nmds_bc_nomis)$sites)
#add columns to data frame 
data.scores$Sample = metadata_nmds$Sample
data.scores$Site = metadata_nmds$Site_c
head(data.scores)

##plot nmds
##stress=0.172395
##k=2
nmds_bc_GFS_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 6, aes(colour = Site))+ 
 #stat_ellipse(aes(x=NMDS1, y=NMDS2,color=Site),type = "norm")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Region", y = "NMDS2", shape = "Type") 

nmds_bc_GFS_plot

##Include Smooth line to show the latitude
##First use
colors<-c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
          "#3EBCB6","#82581FFF","#2F509EFF",
          "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E")

plot(x=data.scores$NMDS1, y=data.scores$NMDS2, type="n", xlim = c(-2, 2), ylim = c(-1.3, 1.3))
points(nmds_bc_nomis, display = "sites", cex = 2.3, pch=19,col=alpha(colors[factor(data.scores$Site)], 0.8))
ordisurf(nmds_bc_nomis, metadata_nmds$lat_sp, add = TRUE, col="blue",labcex=0)

## Statistical analyses
vegan_matrix_GFS<- vegan_otu(prune_Uganda)
GFS_bray <- vegdist(log1p(vegan_matrix_GFS), method="bray")
GFS_ado <- adonis2(GFS_bray ~ Site_c, permutations = 999, method = "bray", data=metadata_nmds)
GFS_ado_latitude <- adonis2(GFS_bray ~ lat_attribute, permutations = 999, method = "bray", data=metadata_nmds)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = GFS_bray ~ Site_c, data = sample.data.frame(prune_Uganda), permutations = 999, method = "bray")
# Df SumOfSqs     R2      F Pr(>F)    
# Site_c     9   16.785 0.3098 6.8324  0.001 ***
#   Residual 137   37.396 0.6902                  
# Total    146   54.182 1.0000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##effect of latitude
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = GFS_bray ~ lat_attribute, data = metadata_nmds, permutations = 999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# lat_attribute   1    4.137 0.07636 11.988  0.001 ***
#   Residual      145   50.044 0.92364                  
# Total         146   54.182 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

library(pairwiseAdonis)
pairwise_GFS <- pairwise.adonis(GFS_bray, metadata_nmds$Site_c, p.adjust.m="holm")
# 
# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1          New_Zealand vs Alps  1 3.1782417 11.512354 0.20738364   0.001      0.045   .
# 2         New_Zealand vs Nepal  1 2.5545641  9.243586 0.21375624   0.001      0.045   .
# 3  New_Zealand vs kyrgyzstan  1 3.5985671 13.784926 0.28256527   0.001      0.045   .
# 4         New_Zealand vs Chile  1 1.1926894  3.910894 0.12255672   0.001      0.045   .
# 5        New_Zealand vs Alaska  1 2.5351775  8.888425 0.21219287   0.001      0.045   .
# 6     New_Zealand vs Greenland  1 1.8442631  6.914576 0.22366718   0.001      0.045   .
# 7      New_Zealand vs Caucasus  1 2.9523560 11.268725 0.23345811   0.001      0.045   .
# 8       New_Zealand vs Ecuador  1 2.0279543  7.574851 0.21908558   0.001      0.045   .
# 9        New_Zealand vs Norway  1 2.2002549  8.167314 0.23224162   0.001      0.045   .
# 10               Alps vs Nepal  1 1.2605281  4.548463 0.10210146   0.001      0.045   .
# 11        Alps vs kyrgyzstan  1 2.5149213  9.524608 0.18851425   0.001      0.045   .
# 12               Alps vs Chile  1 1.5444644  5.134075 0.13119193   0.001      0.045   .
# 13              Alps vs Alaska  1 1.0935921  3.841872 0.08967562   0.001      0.045   .
# 14           Alps vs Greenland  1 1.2183025  4.517612 0.13087847   0.001      0.045   .
# 15            Alps vs Caucasus  1 0.7426476  2.805427 0.06124661   0.002      0.045   .
# 16             Alps vs Ecuador  1 2.2170295  8.204277 0.19911227   0.001      0.045   .
# 17              Alps vs Norway  1 1.4959698  5.508000 0.14303520   0.001      0.045   .
# 18       Nepal vs kyrgyzstan  1 1.9682796  7.556326 0.19598148   0.001      0.045   .
# 19              Nepal vs Chile  1 1.3270923  4.259658 0.15073284   0.001      0.045   .
# 20             Nepal vs Alaska  1 1.3170500  4.573929 0.13623454   0.001      0.045   .
# 21          Nepal vs Greenland  1 1.7389208  6.513499 0.24566728   0.001      0.045   .
# 22           Nepal vs Caucasus  1 1.4061189  5.375591 0.14007839   0.001      0.045   .
# 23            Nepal vs Ecuador  1 2.3164838  8.639917 0.27307015   0.001      0.045   .
# 24             Nepal vs Norway  1 2.0650598  7.646071 0.24949597   0.001      0.045   .
# 25       kyrgyzstan vs Chile  1 1.8742748  6.491945 0.20614621   0.001      0.045   .
# 26      kyrgyzstan vs Alaska  1 2.2084607  8.188597 0.21442519   0.001      0.045   .
# 27   kyrgyzstan vs Greenland  1 1.9527239  8.072358 0.27766436   0.001      0.045   .
# 28    kyrgyzstan vs Caucasus  1 2.3052916  9.361608 0.21589624   0.001      0.045   .
# 29     kyrgyzstan vs Ecuador  1 2.6097633 10.603151 0.30642155   0.001      0.045   .
# 30      kyrgyzstan vs Norway  1 2.4757287  9.982133 0.29374651   0.001      0.045   .
# 31             Chile vs Alaska  1 1.2977613  3.983387 0.14762368   0.001      0.045   .
# 32          Chile vs Greenland  1 1.2728115  3.975568 0.22116509   0.001      0.045   .
# 33           Chile vs Caucasus  1 1.6174920  5.617196 0.17221579   0.001      0.045   .
# 34            Chile vs Ecuador  1 1.2920550  4.136985 0.19572258   0.001      0.045   .
# 35             Chile vs Norway  1 1.4376638  4.564312 0.21166046   0.001      0.045   .
# 36         Alaska vs Greenland  1 1.4184767  5.032398 0.20940058   0.001      0.045   .
# 37          Alaska vs Caucasus  1 1.4587434  5.397730 0.14433310   0.001      0.045   .
# 38           Alaska vs Ecuador  1 2.1784548  7.751633 0.26054478   0.001      0.045   .
# 39            Alaska vs Norway  1 1.7669213  6.241593 0.22100713   0.001      0.045   .
# 40       Greenland vs Caucasus  1 1.1765330  4.800465 0.17267571   0.001      0.045   .
# 41        Greenland vs Ecuador  1 1.1565206  4.740337 0.26720671   0.001      0.045   .
# 42         Greenland vs Norway  1 0.3778133  1.526799 0.10510224   0.058      0.058    
# 43         Caucasus vs Ecuador  1 2.1482909  8.640759 0.24943908   0.001      0.045   .
# 44          Caucasus vs Norway  1 1.4521478  5.800169 0.18239428   0.001      0.045   .
# 45           Ecuador vs Norway  1 1.3568216  5.368108 0.25122055   0.001      0.045   .


library(pairwiseAdonis)
pairwise_GFS_lat <- pairwise.adonis(GFS_bray, metadata_nmds$lat_attribute, p.adjust.m="holm")
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1 FALSE vs TRUE  1  4.137412 11.98788 0.07636181   0.001      0.001  **

