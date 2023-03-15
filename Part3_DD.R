##Distance Decay - NOMIS - Biogeography
library(speedyseq)
library(phyloseq)
library(phyloseqCompanion)
library(geosphere)
library(vegan)
library(rbiom)
library(scales)
library(ggplot2)
library(dplyr)
library(fishualize)
library(ggpubr)
library(reshape2)
library(performance)
library(fANCOVA)

### Import all data 
tree_NOMIS<- read_tree("20220411_NOMIS_TREE_madebyG_rarefied.tree")
NOMIS_FR <- readRDS("NOMIS_FR_20221026.RData")
tax<-read.csv(file="1022_NOMIS_16S_merged_taxonomy.csv",sep=",",header=TRUE,row.names=1)
tax_NOMIS <- tax_table(as.matrix(tax))
otu_nomis_full <- otu_table(NOMIS_FR, taxa_are_rows=T)
metadata_glaciers="metadata_glaciers.tsv"
metadata_glaciers<-import_qiime_sample_data(metadata_glaciers)

##Merge everything
merge_nomis_data <- merge_phyloseq(otu_nomis_full, tax_NOMIS, metadata_glaciers, tree_NOMIS)
##Prune Uganda
uganda="Uganda"
prune_Uganda <- subset_samples(merge_nomis_data, !Site_c %in% uganda)
nomis_metadata <- sample.data.frame(prune_Uganda)

vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

############################################################Distance-Decay patterns on Bray-Curtis
##Full distance decay with all samples
##building dataframe encompassing latitude and longitude from all glaciers
NOMIS_df <- data.frame(nomis_metadata$lon_sp, nomis_metadata$lat_sp)
NOMIS_df$nomis_metadata.lon_sp<-as.numeric(NOMIS_df$nomis_metadata.lon_sp)
NOMIS_df$nomis_metadata.lat_sp<-as.numeric(NOMIS_df$nomis_metadata.lat_sp)

dist_geo_all<-distm(NOMIS_df, NOMIS_df, fun=distGeo)
dist_geo_all <- as.matrix(dist_geo_all)
diag(dist_geo_all)=NA
dist_geo_all_diag<-t(matrix(t(dist_geo_all)[which(!is.na(dist_geo_all))],nrow=146,ncol=147))
#write.csv(dist_geo_all_diag, "dist_geo_nomis.csv")

vegan_matrix_all <- vegan_otu(prune_Uganda)
allregion_bray <- vegdist(log1p(vegan_matrix_all), method="bray")
allregion_m<- as.matrix(allregion_bray)
diag(allregion_m)=NA
allregion_diag<-t(matrix(t(allregion_m)[which(!is.na(allregion_m))],nrow=146,ncol=147))
min(allregion_diag)

mantel(allregion_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)

# Mantel statistic based on Pearson's product-moment correlation 
# # 
# Call:
# mantel(xdis = allregion_diag, ydis = dist_geo_all_diag, method = "pearson",      permutations = 999, na.rm = TRUE) 
# 
# Mantel statistic r: 0.5306 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0335 0.0458 0.0563 0.0705 
# Permutation: free
# Number of permutations: 999

dist.all.bray <- data.frame(BC.dist_bc=as.vector(allregion_diag), BC.sim_bc=as.vector(1-allregion_diag), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers")
ggplot(dist.all.bray, aes(x=geo.dist/1000, y=BC.sim_bc)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
 # facet_wrap(~Method) +
  stat_smooth(method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Bray-Curtis") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()

##Getting the slope coefficients
  dist.all_bray %>% 
  do({
    mod = lm(BC.sim_bc ~ geo.dist, data = dist.all_bray)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
  
# Intercept         Slope
#   (Intercept) 0.2198439 -9.078998e-09

#################################################Full distance decay using Weighted-Unifrac distance
##Full distance decay with all samples
##building dataframe encompassing latitude and longitude from all glaciers
asv_table_unif <-otu_table(prune_Uganda, taxa_are_rows=T)
unifrac_nomis <- rbiom::unifrac(asv_table_unif, weighted=T, tree=phy_tree(prune_Uganda))
allregion_unif<- as.matrix(unifrac_nomis)
diag(allregion_unif)=NA
allregion_diag_unif<-t(matrix(t(allregion_unif)[which(!is.na(allregion_unif))],nrow=146,ncol=147))

mantel(allregion_diag_unif, dist_geo_all_diag, method = "pearson", permutations = 9999, na.rm = TRUE)

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = allregion_diag_unif, ydis = dist_geo_all_diag,      method = "pearson", permutations = 9999, na.rm = TRUE) 
# 
# Mantel statistic r: 0.1389 
#       Significance: 5e-04 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0455 0.0607 0.0724 0.0881 
# Permutation: free
# Number of permutations: 9999

dist.all.wunif <- data.frame(dis_wunif=as.vector(allregion_diag_unif), sim_wunif=as.vector(1-allregion_diag_unif), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers")

ggplot(dist.all.wunif, aes(x=geo.dist/1000, y=sim_wunif)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
 # facet_wrap(~Method) +
  stat_smooth(method="lm", formula= y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Weighted Unifrac") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()

dist.all.wunif %>% 
  do({
    mod = lm(sim_wunif ~ geo.dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
# Intercept         Slope
# (Intercept) 0.6665992 -2.293525e-09

####################################################Full distance decay using Unweighted unifrac distance
asv_table_unif <-otu_table(prune_Uganda, taxa_are_rows=T)
uwnifrac_nomis <- rbiom::unifrac(asv_table_unif, weighted=F, tree=phy_tree(prune_Uganda))
allregion_uwnif<- as.matrix(uwnifrac_nomis)
diag(allregion_uwnif)=NA
allregion_diag_uwnif<-t(matrix(t(allregion_uwnif)[which(!is.na(allregion_uwnif))],nrow=146,ncol=147))

mantel(allregion_diag_uwnif, dist_geo_all_diag, method = "pearson", permutations = 9999, na.rm = TRUE)

# 
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = allregion_diag_uwnif, ydis = dist_geo_all_diag,      method = "pearson", permutations = 9999, na.rm = TRUE) 
# 
# Mantel statistic r: 0.2841 
#       Significance: 1e-04 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0409 0.0533 0.0646 0.0789 
# Permutation: free
# Number of permutations: 9999

dist.all.uw <- data.frame(dis_uw=as.vector(allregion_diag_uwnif), sim_uw=as.vector(1-allregion_diag_uwnif), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers")

ggplot(dist.all.uw, aes(x=geo.dist/1000, y=sim_uw)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  # facet_wrap(~Method) +
  stat_smooth(method="lm", formula= y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - UnWeighted Unifrac") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()


dist.all.uw %>% 
  do({
    mod = lm(sim_uw ~ geo.dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })
# Intercept         Slope
# (Intercept) 0.3534007 -3.425515e-09

##################################################### Full distance decay using SORENSEN
vegan_matrix_all <- vegan_otu(prune_Uganda)
allregion_sor <- vegdist(vegan_matrix_all, binary=T)
allregion_sor_m<- as.matrix(allregion_sor)
diag(allregion_sor_m)=NA
allregion_sor_diag<-t(matrix(t(allregion_sor_m)[which(!is.na(allregion_sor_m))],nrow=146,ncol=147))
min(allregion_sor_diag)

mantel(allregion_sor_diag, dist_geo_all_diag, method = "pearson", permutations = 999, na.rm = TRUE)
# 
# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = allregion_sor_diag, ydis = dist_geo_all_diag, method = "pearson",      permutations = 999, na.rm = TRUE) 
# 
# Mantel statistic r: 0.5492 
#       Significance: 0.001 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0351 0.0462 0.0537 0.0734 
# Permutation: free
# Number of permutations: 999

dist.all.sor <- data.frame(dis_sor=as.vector(allregion_sor_diag), sim_sor=as.vector(1-allregion_sor_diag), geo.dist=as.vector(dist_geo_all_diag), Method="All glaciers")
ggplot(dist.all.sor, aes(x=geo.dist/1000, y=sim_sor)) +
  geom_point(size=1.2) + ylim(0,1) + xlim (0,20000)+
  stat_smooth(method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  ylab("Community Similarity - Sorensen") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  theme(strip.text.x=element_text(size=10, color="black", face="bold.italic")) +
  theme(strip.background=element_rect(colour="black", fill="white")) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        legend.text=element_text(size=14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_log10()

dist.all.sor %>% 
  do({
    mod = lm(sim_sor ~ geo.dist, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  })

# Intercept         Slope
# (Intercept) 0.2240334 -9.032569e-09


##Combine different decay patterns together
##Lets combine BC and Sorensen
melt_bc <- melt(allregion_diag)
melt_bc$dis <- "BC"
melt_sor <- melt(allregion_sor_diag)
melt_sor$dis <- "SOR"
melt_geo <- melt(dist_geo_all_diag)
colnames(melt_geo) <- c("Var1","Var2","dist_geo")
bind_bcsor <- rbind(melt_bc, melt_sor)
merge_dis_geo <- merge(bind_bcsor, melt_geo)

##Plot everything!
bcsor_plot<- ggplot(merge_dis_geo) + 
  geom_point(aes(y=1-value, x = dist_geo, color= dis)) + 
  stat_smooth(aes(y=1-value, x = dist_geo, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_fish(option="Scarus_quoyi",discrete=T)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
bcsor_plot


##Lets combine unweighted and weighted unifrac
melt_unifrac <- melt(allregion_diag_unif)
melt_unifrac$dis <- "wunifrac"
melt_uwunifrac <- melt(allregion_diag_uwnif)
melt_uwunifrac$dis <- "unwunifrac"
bind_unifrac<- rbind(melt_unifrac, melt_uwunifrac)
merge_dis_unif <- merge(bind_unifrac, melt_geo)

##Plot everything!
unifrac_plot<- ggplot(merge_dis_unif) + 
  geom_point(aes(y=1-value, x = dist_geo, color= dis)) + 
  stat_smooth(aes(y=1-value, x = dist_geo, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  scale_color_fish(option="Trimma_lantana",discrete=T)+
  ylab("Community Similarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
unifrac_plot


##Create a plot with all the distance matrices
bind_all_dist <- rbind(melt_unifrac,melt_uwunifrac,melt_bc,melt_sor)
merge_alldist_geo <- merge(bind_all_dist, melt_geo)

get_fishcol <- fish(4, option="Trimma_lantana")

alldist_plot<- ggplot(merge_alldist_geo) + 
  geom_point(aes(y=value, x = dist_geo/1000, color= dis), alpha=0.1) + 
  stat_smooth(aes(y=value, x = dist_geo/1000, color=dis), method="lm", formula=y ~ (x), size=1.2, se=FALSE, linetype="solid") +
  #scale_color_manual(values=c("#F74A24FF","#4AB2B8FF","#F5AE21FF","#791214FF"))+
  scale_color_fish(option="Scarus_quoyi",discrete=T,direction=-1)+
  stat_regline_equation(aes(y=value, x = dist_geo/1000, color=dis, label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y ~ (x)
  ) +
  ylab("Community Dissimilarity") +
  xlab("Geographic distance (km)") + ggtitle("") + theme_bw() +
  scale_x_log10() +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggsave(alldist_plot,path='dissimilarity.pdf', device='pdf', dpi=100)
  dev.off
  
  # save the plot in the temporary folder
  file_temp <- tempfile(fileext = '.pdf')
  ggsave(file.path, plot = alldist_plot, "dissimilarity2023.pdf",dpi=100, width=5, height=7, device="pdf")
  file.path = '/Users/ezzat/Google Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_October2022'
  folder_perc ="/Users/ezzat/Google Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_October2022"
  file.copy(file_temp, file.path(folder_perc, "dissimilarity.pdf"))
  
###Distance-Decay patterns comparisons
##Bray-Curtis vs Sorensen
##Lets try another way
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)

  ## the covariate -> x
  ggscatter(
   manco, x = "geo.dist.x", y = "value",
    color = "variable", add = "reg.line"
  )+
    stat_regline_equation(
      aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = variable)
    )
 ## there is a linear relationship! 

## Homogeneity of regression slopes  
#This assumption checks that there is no significant interaction between the covariate and the grouping variable. 
#This can be evaluated as follow:
  
merge_bcsor <- merge(dist.all.bray, dist.all.sor, by="row.names")
anco_bcsor <- merge_bcsor[c("BC.dist_bc","dis_sor","geo.dist.x")]
manco <- melt(anco_bcsor,  id=c("geo.dist.x")) 
  
##Weighted vs Unweighted unifrac
merge_unifrac <- merge(dist.all.wunif, dist.all.uw, by="row.names")
anco_unifrac <- merge_unifrac[c("dis_wunif","dis_uw","geo.dist.x")]
manco_unifrac<- melt(anco_unifrac, id=c("geo.dist.x")) 
  
    
manco %>% anova_test(sqrt(value) ~ variable*geo.dist.log)  
# ANOVA Table (type II tests)
# 
# Effect DFn   DFd         F        p p<.05      ges
# 1            variable   1 42920    36.292 1.71e-09     * 8.45e-04
# 2          geo.dist.x   1 42920 17631.160 0.00e+00     * 2.91e-01
# 3 variable:geo.dist.x   1 42920     0.116 7.34e-01       2.70e-06 

## there is homogeneity of regression slopes as the interaction term was not significant! p=0.73 

##Normality of residuals
# Fit the model, the covariate goes first
model <- lm((value) ~ geo.dist.log + variable, data = manco)
# Inspect the model diagnostic metrics
model.metrics <- augment(model) %>%
  select(-.hat, -.sigma, -.fitted) # Remove details
head(model.metrics, 3)
# Assess normality of residuals using shapiro wilk test
## our sample size is too big, we have to use anaother test.
#Anderson-Darling normality test
library(nortest)
ad.test(model$residuals)
# Anderson-Darling normality test
# 
# data:  model.metrics$.resid
# A = 196.93, p-value < 2.2e-16

##Homogeneity of regression slopes
manco_unifrac %>% anova_test(log10(value) ~variable*geo.dist.x)
# ANOVA Table (type II tests)
# Effect DFn   DFd          F        p p<.05   ges
# 1            variable   1 32000 120225.886 0.00e+00     * 0.790
# 2          geo.dist.x   1 32000   1604.933 0.00e+00     * 0.048
# 3 variable:geo.dist.x   1 32000    436.514 2.73e-96     * 0.013
##Violation of the assumption, we try with non-parametric ANCOVA


library(effectsize)
manco$geo.dist.log <- log((manco$geo.dist.x)/1000)
glm_dd_bcsor <- (lm(manco, formula = value ~ geo.dist.log + geo.dist.log:variable + variable))
check_model(glm_dd_bcsor)
# 
# Call:
#   lm(formula = value ~ geo.dist.log + geo.dist.log:variable + variable, 
#      data = manco)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.53538 -0.04008  0.00924  0.05164  0.27154 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   0.5615632  0.0026189 214.423   <2e-16 ***
#   geo.dist.log                  0.0348620  0.0003041 114.655   <2e-16 ***
#   variabledis_sor              -0.0005688  0.0037038  -0.154    0.878    
# geo.dist.log:variabledis_sor -0.0004744  0.0004300  -1.103    0.270    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07383 on 42920 degrees of freedom
# Multiple R-squared:  0.377,	Adjusted R-squared:  0.377 
# F-statistic:  8659 on 3 and 42920 DF,  p-value: < 2.2e-16

manco_unifrac$geo.dist.log <- log((manco_unifrac$geo.dist.x)/1000)
glm_dd_wufrac <- lm(manco_unifrac, formula = value ~ geo.dist.log + geo.dist.log:variable + variable)
check_model(glm_dd_wufrac)

# Call:
#   lm(formula = value ~ geo.dist.log + geo.dist.log:variable + variable, 
#      data = manco_unifrac)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.38057 -0.05302 -0.00050  0.04985  0.32084 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 0.2377637  0.0027120  87.670  < 2e-16 ***
#   geo.dist.log                0.0135885  0.0003149  43.156  < 2e-16 ***
#   variabledis_uw              0.2945001  0.0038354  76.785  < 2e-16 ***
#   geo.dist.log:variabledis_uw 0.0033344  0.0004453   7.488 7.12e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07645 on 42920 degrees of freedom
# Multiple R-squared:  0.8203,	Adjusted R-squared:  0.8203 
# F-statistic: 6.53e+04 on 3 and 42920 DF,  p-value: < 2.2e-16


