##Computing between and within BC and Sorensen indices of diversity from data that were filtered, grouped within the same patch and rarefied 
library(vegan)
library(adespatial)
library(fishualize)
merge_glaciers_data <- as.data.frame(sample_data(prune_Uganda))
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 45737 taxa and 147 samples ]:
#   sample_data() Sample Data:        [ 147 samples by 38 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 45737 taxa by 7 taxonomic ranks ]:
#   taxa are rows

vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}
vegan_matrix_nomis <- vegan_otu(prune_Uganda)
##logtransform
vegan_matrix_nomis <- (vegan_matrix_nomis)
asv_region_bray<- vegdist(log1p(vegan_matrix_nomis), method="bray") ## Bray-Curtis ## log transformation downweight the importance of abundant species.
asv_region_sorensen <-vegdist(vegan_matrix_nomis, binary=T) ## Sorensen
asv_region_jaccard<- vegdist((vegan_matrix_nomis), method="jaccard") ## Bray-Curtis ## log transformation downweight the importance of abundant species.

##meandist
nomis_meandist_bray <- vegan::meandist(asv_region_bray, merge_glaciers_data$Site_c)
hist(meandist_bray_within)
meandist_bray_within <- c(0.764256,0.735157,0.696905,0.846161,0.699021,0.679834,0.680074,0.73785,0.728873,0.705885)
nomis_meandist_sorensen <- vegan::meandist(asv_region_sorensen, merge_glaciers_data$Site_c)

##sddist
sddist <- function (dist, grouping, ...) 
{
  mergenames <- function(X, Y, ...) {
    xy <- cbind(X, Y)
    xy <- apply(xy, 1, sort)
    apply(xy, 2, paste, collapse = " ")
  }
  grouping <- factor(grouping, exclude = NULL)
  cl <- outer(grouping, grouping, mergenames)
  cl <- cl[lower.tri(cl)]
  n <- table(grouping)
  take <- matrix(TRUE, nlevels(grouping), nlevels(grouping))
  diag(take) <- n > 1
  take[upper.tri(take)] <- FALSE
  out <- matrix(NA, nlevels(grouping), nlevels(grouping))
  out[take] <- tapply(dist, cl, sd)
  out[upper.tri(out)] <- t(out)[upper.tri(out)]
  rownames(out) <- colnames(out) <- levels(grouping)
  class(out) <- c("meandist", "matrix")
  attr(out, "n") <- table(grouping)
  out
}

nomis_sddist_bray <- with(merge_glaciers_data, sddist(asv_region_bray, Site_c))
nomis_sddis_sorensen <- with(merge_glaciers_data, sddist(asv_region_sorensen, Site_c))
##Plotting Dissimilarity
library(RColorBrewer)
nomis_braycurtis <- as.matrix(nomis_meandist_bray)
nomis_sorensen <- as.matrix(nomis_meandist_sorensen)

# Heatmap
fish_color<- fish(n=5, option = "Zebrasoma_velifer", end=1,  begin=0.3)
fish_color_bray<- fish(n=5, option = "Coris_gaimard",direction=1,begin=0.2)
heatmap(nomis_braycurtis, Rowv=NA, Colv=NA, col = fish_color)
legend(x="bottomright", legend=c("min", "ave", "max"), 
       col=c("#008AC2","#3F459B"))
library(pheatmap)

##Sorensen 
pheatmap(nomis_sorensen,
         display_numbers = T,
         color =fish_color, 
         fontsize_number = 8,
         treeheight_row=0,
         treeheight_col=0, cluster_rows=F, cluster_cols=F, na_col="white")

##Bray-Curtis 
pheatmap(nomis_braycurtis,
         display_numbers = T,
         color =fish_color_bray, 
         fontsize_number = 8,
         treeheight_row=0,
         treeheight_col=0,cluster_rows=F, cluster_cols=F, na_col="white")

nomis_braycurtis[lower.tri(nomis_braycurtis)] <- NA
pheatmap(nomis_braycurtis, cluster_rows=F, cluster_cols=F, na_col="white")
##Plotting the boxplot of BC using mean and sd

# box_bc <- ggplot(boxplot_data,aes(x=Site))+geom_boxplot(aes(lower=Mean_BC-SD,upper=Mean_BC+SD,middle=Mean_BC,ymin=Mean_BC-3*SD,ymax=Mean_BC+3*SD),stat="identity")
# box_bc<- box_bc+ scale_y_continuous(breaks = seq(0, 1, by = 0.2))+ theme_bw()
# 

asv_region_bray_mat <- as.matrix(asv_region_bray)
asv_region_sor_mat <- as.matrix(asv_region_sorensen)
dissimilarity <- data.frame(b=as.vector(asv_region_bray_mat[upper.tri(asv_region_bray_mat)]),
                  s=as.vector(asv_region_sor_mat[upper.tri(asv_region_sor_mat)]))
ggplot(dissimilarity, aes(x=s, y=b))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

fit_dissimilarity=lm(b~s,data=dissimilarity)
summary(fit_dissimilarity)

# Call:
#   lm(formula = b ~ s, data = dissimilarity)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.078282 -0.005554 -0.000591  0.005519  0.053011 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.0245362  0.0009484  -25.87   <2e-16 ***
#   s            1.0341863  0.0011072  934.06   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01051 on 10729 degrees of freedom
# Multiple R-squared:  0.9879,	Adjusted R-squared:  0.9879 
# F-statistic: 8.725e+05 on 1 and 10729 DF,  p-value: < 2.2e-16

