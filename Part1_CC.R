### Comparisons with other samples from Cryospheric Ecosystems from Massimo's paper

library(phyloseq)
library(data.table)
library(tidyverse)
library(Rmisc)
library(vegan)
library(wesanderson)
library(paletteer)
library(FSA)

## read in otu table 
asv_table_NOMISCRYO = fread("PP1_table.tsv", header=T)

## remove columns containing GFSs from previous expeditions 
Data_cryo <- subset(asv_table_NOMISCRYO, select = -c(GL1_UP_1A_16S_sed, GL17_UP_1A_16S_sed, GL19_UP_1A_16S_sed, GL20_UP_1A_16S_sed, GL43_UP_1_16S_sed, GL51_UP_3_16S_sed, GL54_UP_1_16S_sed, GL56_UP_1_16S_sed, GL6_UP_1A_16S_sed, GL7_UP_1A_16S_sed))

## Load metadata CRYO
metadata_cryo="/PP1_metadata_date.tsv"
metadata_cryo<-import_qiime_sample_data(metadata_cryo)
## Remove the taxonomy 
Data_cryo <- subset(Data_cryo, select = -c(Taxonomy))

## We start here with an ASV table that takes into account the 148 samples and the ASVs that were filtered -- We keep samples from Uganda
merge_nomis_cryo <- merge(NOMIS_full_asv_sedimentup, Data_cryo, by.x="row.names", by.y="ASV", all.y=T, all.x=T)
## Now we need to keep only data from Cryospheric ecosystems
## Load metadata CRYO + NOMIS
metadata_merge_cryonom="metadata_merge_nomis_cryo.txt"
metadata_cryo_nomis <-import_qiime_sample_data(metadata_merge_cryonom)
##Replace NAs by 0 in merge_nomis_cryo because those NA are taxa that were absent from the GFS
merge_nomis_cryo[is.na(merge_nomis_cryo)] <- 0

#saveRDS(merge_nomis_cryo,"merge_nomis_cryo.RDS")
merge_nomis_cryo <- readRDS("merge_nomis_cryo.RDS")

asv_table_NOMISCRYO = as.data.frame(merge_nomis_cryo)
unAsIs <- function(X) {
   if("AsIs" %in% class(X)) {
     class(X) <- class(X)[-match("AsIs", class(X))]
   }
  X
}

##Remove AsIs --- tell that the first column of ASVs are the rownames
asv_table_NOMISCRYO$Row.names<-unAsIs(asv_table_NOMISCRYO$Row.names)
rownames(asv_table_NOMISCRYO) <- asv_table_NOMISCRYO$Row.names
asv_table_NOMISCRYO$Row.names <- NULL
asv_table_NOMISCRYO_m <- as.matrix(asv_table_NOMISCRYO)
## transform as otu table --> yeahhh
asv_tableNOMISCRYO_final <- otu_table(asv_table_NOMISCRYO_m, taxa_are_rows=T)

##now merge asv_table and metadata
merge_asvnomiscryo_metadata <- merge_phyloseq(asv_tableNOMISCRYO_final, metadata_cryo_nomis)

##Remove everything else that is not cryospheric!!!
eco = c("Yes")
Cryosphere_sample <- subset_samples(merge_asvnomiscryo_metadata, Cryosphere %in% eco)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 559355 taxa and 433 samples ]
# sample_data() Sample Data:       [ 433 samples by 9 sample variables ]
##Ouhouuuu ça marche!!!

library(R.filesets)
#saveRDS(Cryosphere_sample, "Cryosphere_sample.RDS")

##We need to compute a rarefaction curve for these samples..because the may not all fit well.
prune_nomiscryo <- prune_samples(sample_sums(Cryosphere_sample) >10000,Cryosphere_sample)
NOMIS_CRYOS<- rarefy_even_depth(prune_nomiscryo, sample.size=min(sample_sums(prune_nomiscryo)), rngseed=678, replace=F, trimOTUs=TRUE, verbose=TRUE)

#saveRDS(NOMIS_CRYOS,file="20232402_NOMISCRYOS.RData")
NOMIS_CRYOS<-readRDS("20232402_NOMISCRYOS.RData")
#load("NOMISCROYS.RData")
meta <- sample.data.frame(NOMIS_CRYOS)
##Dans la table à la base j'ai remplacé ce qui est permafrost par "terrestrial"
NOMIS_CRYOS_sub_rock <- subset_samples(NOMIS_CRYOS, Habitat != "Rock")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub_rock, Habitat != "Ice/Snow")
NOMIS_CRYOS_sub <- subset_samples(NOMIS_CRYOS_sub, Habitat !="Water_lake")

##Calculate ASV richness
diversity_cryosphere <- plot_richness(NOMIS_CRYOS_sub, measures=c("Observed","Shannon"), color="Habitat", shape="Ecosystem")
alphadt_cryosphere <-data.table(diversity_cryosphere$data)

OR_cryo <- subset(alphadt_cryosphere, variable == "Observed")
ASVrichness_cryo <- OR_cryo %>% 
  group_by(Habitat) %>% 
  summarise(average=mean(value), std=sd(value))

cryosphere_richness <- ggplot(OR_cryo,aes(x=Habitat,y=value, color=Habitat)) + 
  geom_boxplot() + 
  geom_jitter(aes(group=Habitat), position=position_jitterdodge()) + 
  scale_fill_viridis(discrete = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

median_cryosphere <- alphadt_cryosphere %>%
  group_by(Habitat) %>% summarize(median=median(value),x = quantile(value, c(0.25, 0.5, 0.75)))

# # A tibble: 27 × 3
# # Groups:   Habitat [9]
# Habitat       median     x
# <chr>          <dbl> <dbl>
# 1 Cryoconite       126   93 
# 2 Cryoconite       126  126 
# 3 Cryoconite       126  405 
# 4 GFS             1267  987 
# 5 GFS             1267 1267 
# 6 GFS             1267 1551.
# 7 Glacier          542  468.
# 8 Glacier          542  542 
# 9 Glacier          542  678.
# 10 Microbial mat    293  259 
# # … with 17 more rows

Shannon_cryo <- subset(alphadt_cryosphere, variable == "Shannon")
Shannondiv_cryo <- Shannon_cryo %>% 
  group_by(Habitat) %>% 
  summarize(median=median(value),x = quantile(value, c(0.25, 0.5, 0.75)))
# 
# # A tibble: 27 × 3
# # Groups:   Habitat [9]
# Habitat       median     x
# <chr>          <dbl> <dbl>
# 1 Cryoconite      3.39  3.08
# 2 Cryoconite      3.39  3.39
# 3 Cryoconite      3.39  4.59
# 4 GFS             6.01  5.50
# 5 GFS             6.01  6.01
# 6 GFS             6.01  6.31
# 7 Glacier         4.95  4.81
# 8 Glacier         4.95  4.95
# 9 Glacier         4.95  5.34
# 10 Microbial mat   4.75  4.27

cryosphere_shannon <- ggplot(Shannon_cryo,aes(x=Habitat,y=value, color=Habitat)) + 
  geom_boxplot()+ 
  geom_jitter(aes(group=Habitat), position=position_jitterdodge()) + 
  scale_fill_viridis(discrete = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# ##Lets try a PCOA
# library(phyloseqCompanion)
# metadata_pcoa_cryo <- sample.data.frame(sample_data(NOMIS_CRYOS_sub))
# asv_table_cryo <- as.matrix((log1p(otu_table(Cryosphere, taxa_are_rows=T))))
# 
# ##Compute PCoA
# nomis_bray_dist_cryo <- matrix(vegdist((asv_table_cryo), method = "bray"))
# pcoa_cryo<-ecodist::pco(nomis_bray_dist_cryo)
# 
# pcoa_bc_cryo_df <- data.frame(pcoa1 = pcoa_cryo$vectors[,1], 
#                                pcoa2 = pcoa_cryo$vectors[,2])
# merge_data_pcoa_cryo <- cbind(pcoa_bc_cryo_df,
#                          Site = metadata_pcoa_cryo$Habitat)
# ##PCoA plot
# pcoa_bc_plot <- ggplot(data = merge_data_pcoa_cryo, aes(x=pcoa1, y=pcoa2, color=Site)) +
#   geom_point() +
#   labs(x = "PC1",
#        y = "PC2", 
#        title = "PCoA Bray-Curtis") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))  
# 
# pcoa_bc_plot

metadata_nmds_cryo <- sample.data.frame(sample_data(NOMIS_CRYOS_sub))
ord <- ordinate(NOMIS_CRYOS_sub, method = "NMDS", distance = "bray", trymax = 999, k=2)
stressplot(ord)

data.scores = as.data.frame(scores(ord)$sites)

#add columns to data frame 
data.scores$Sample = metadata_nmds_cryo$Sample
data.scores$Site = metadata_nmds_cryo$Habitat
head(data.scores)

##plot NMDS
##stress=0.1103596
##k=2
nmds_bc_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 6, aes(color = Site), alpha=0.6) + 
  #stat_ellipse(aes(x=NMDS1, y=NMDS2,color=Site),type = "norm")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Biomes", y = "NMDS2", shape = "Type") 

## Statistical analyses to compare diversity across different biomes
library(performance)
library(multcomp)
habitats_richness <-glm(sqrt(value) ~ as.factor(Habitat), data = OR_cryo)
summary(habitats_richness)
aov_region <-anova(habitats_richness)
check_model(habitats_richness)
check_normality(habitats_richness)
check_heteroscedasticity(habitats_richness)
##ollowing normality // Homoscedasticity of variances
habitats_richness_np <- aov(value ~ Habitat, data = OR_cryo)
summary(habitats_richness_np)
# > summary(habitats_richness_np)
# Df    Sum Sq  Mean Sq F value Pr(>F)    
# Habitat       8 166669962 20833745   195.4 <2e-16 ***
#   Residuals   399  42547248   106635                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

OR_tukey <- TukeyHSD(habitats_richness_np, "Habitat",ordered=T, method="holm")
# 
# $Habitat
# diff         lwr       upr     p adj
# Water-Permafrost           104.44444  -86.211945  295.1008 0.7410706
# Cryoconite-Permafrost      122.17417 -103.848585  348.1969 0.7549759
# Microbial mat-Permafrost   250.75214  -69.940087  571.4444 0.2660050
# Snow-Permafrost            267.64444   38.105130  497.1838 0.0094004
# Glacier-Permafrost         411.89899   69.333366  754.4646 0.0062554
# Sediment-Permafrost        927.44444  484.799436 1370.0895 0.0000000
# GFS-Permafrost            1103.81064  929.580721 1278.0406 0.0000000
# Terrestrial-Permafrost    2019.05420 1799.165250 2238.9432 0.0000000
# Cryoconite-Water            17.72973 -185.577640  221.0371 0.9999991
# Microbial mat-Water        146.30769 -158.800651  451.4160 0.8572833
# Snow-Water                 163.20000  -44.009787  370.4098 0.2567855
# Glacier-Water              307.45455  -20.568023  635.4771 0.0862731
# Sediment-Water             823.00000  391.511642 1254.4884 0.0000002
# GFS-Water                  999.36620  855.826432 1142.9060 0.0000000
# Terrestrial-Water         1914.60976 1718.144112 2111.0754 0.0000000
# Microbial mat-Cryoconite   128.57796 -199.793063  456.9490 0.9515409
# Snow-Cryoconite            145.47027  -94.680359  385.6209 0.6215935
# Glacier-Cryoconite         289.72482  -60.039726  639.4894 0.1968362
# Sediment-Cryoconite        805.27027  357.030790 1253.5098 0.0000014
# GFS-Cryoconite             981.63647  793.646857 1169.6261 0.0000000
# Terrestrial-Cryoconite    1896.88003 1665.935924 2127.8241 0.0000000
# Snow-Microbial mat          16.89231 -313.909052  347.6937 1.0000000
# Glacier-Microbial mat      161.14685 -256.097020  578.3907 0.9552551
# Sediment-Microbial mat     676.69231  174.023583 1179.3610 0.0010957
# GFS-Microbial mat          853.05850  557.935974 1148.1810 0.0000000
# Terrestrial-Microbial mat 1768.30206 1444.122494 2092.4816 0.0000000
# Glacier-Snow               144.25455 -207.792673  496.3018 0.9372719
# Sediment-Snow              659.80000  209.777066 1109.8229 0.0002203
# GFS-Snow                   836.16620  643.962913 1028.3695 0.0000000
# Terrestrial-Snow          1751.40976 1517.022923 1985.7966 0.0000000
# Sediment-Glacier           515.54545   -1.352554 1032.4435 0.0512017
# GFS-Glacier                691.91165  373.156240 1010.6671 0.0000000
# Terrestrial-Glacier       1607.15521 1261.322741 1952.9877 0.0000000
# GFS-Sediment               176.36620 -248.119856  600.8523 0.9322450
# Terrestrial-Sediment      1091.60976  646.431708 1536.7878 0.0000000
# Terrestrial-GFS            915.24356  734.675134 1095.8120 0.0000000

##Shannon diversity
habitats_shannon <-glm(log(value) ~ as.factor(Habitat), data = Shannon_cryo)
summary(habitats_shannon)
aov_region <-anova(habitats_shannon)
check_model(habitats_shannon)
check_normality(habitats_shannon)
check_heteroscedasticity(habitats_shannon)
## Distribution is not normal and problem with homoscedasticity of variances
## using Kruskal-Wallis

habitats_shannon_np <- kruskal.test(Shannon_cryo$value, Shannon_cryo$Habitat)
habitats_shannon_np
# Kruskal-Wallis rank sum test
# 
# data:  Shannon_cryo$value and Shannon_cryo$Habitat
# Kruskal-Wallis chi-squared = 286.02, df = 8, p-value < 2.2e-16

shannon_dunn <- dunnTest(value ~ Habitat,
                    data=Shannon_cryo,
                    method="holm")

# Comparison           Z      P.unadj        P.adj
# 1             Cryoconite - GFS  -8.8178248 1.167046e-18 3.617844e-17
# 2         Cryoconite - Glacier  -2.5387168 1.112598e-02 2.113937e-01
# 3                GFS - Glacier   2.4147240 1.574714e-02 2.677014e-01
# 4   Cryoconite - Microbial mat  -1.8928042 5.838392e-02 8.173748e-01
# 5          GFS - Microbial mat   3.5108041 4.467534e-04 9.381822e-03
# 6      Glacier - Microbial mat   0.6385021 5.231469e-01 1.000000e+00
# 7      Cryoconite - Permafrost  -1.3038648 1.922797e-01 1.000000e+00
# 8             GFS - Permafrost   7.8227455 5.168354e-15 1.498823e-13
# 9         Glacier - Permafrost   1.7317850 8.331186e-02 9.997423e-01
# 10  Microbial mat - Permafrost   1.0191670 3.081237e-01 1.000000e+00
# 11       Cryoconite - Sediment  -4.7395432 2.142005e-06 5.783414e-05
# 12              GFS - Sediment  -1.0996615 2.714797e-01 1.000000e+00
# 13          Glacier - Sediment  -2.3921494 1.675002e-02 2.680003e-01
# 14    Microbial mat - Sediment  -2.9898585 2.791067e-03 5.582134e-02
# 15       Permafrost - Sediment  -4.1336675 3.570199e-05 8.211458e-04
# 16           Cryoconite - Snow  -0.1555476 8.763896e-01 8.763896e-01
# 17                  GFS - Snow   8.4301608 3.451905e-17 1.035572e-15
# 18              Glacier - Snow   2.4161482 1.568567e-02 2.823421e-01
# 19        Microbial mat - Snow   1.7659758 7.739991e-02 1.000000e+00
# 20           Permafrost - Snow   1.1211512 2.622235e-01 1.000000e+00
# 21             Sediment - Snow   4.6377537 3.522161e-06 9.157618e-05
# 22    Cryoconite - Terrestrial -10.4022429 2.421648e-25 8.475769e-24
# 23           GFS - Terrestrial  -4.1240721 3.722323e-05 8.189110e-04
# 24       Glacier - Terrestrial  -4.3789513 1.192518e-05 2.862043e-04
# 25 Microbial mat - Terrestrial  -5.4932352 3.946367e-08 1.104983e-06
# 26    Permafrost - Terrestrial  -9.5849907 9.246601e-22 2.958912e-20
# 27      Sediment - Terrestrial  -0.6242137 5.324872e-01 1.000000e+00
# 28          Snow - Terrestrial -10.0900796 6.111994e-24 2.016958e-22
# 29          Cryoconite - Water  -0.8317863 4.055296e-01 1.000000e+00
# 30                 GFS - Water  10.3703052 3.384422e-25 1.150703e-23
# 31             Glacier - Water   2.1914493 2.841929e-02 4.262894e-01
# 32       Microbial mat - Water   1.4828627 1.381109e-01 1.000000e+00
# 33          Permafrost - Water   0.6587497 5.100565e-01 1.000000e+00
# 34            Sediment - Water   4.5316220 5.853252e-06 1.463313e-04
# 35                Snow - Water  -0.6358456 5.248771e-01 1.000000e+00
# 36         Terrestrial - Water  11.3670173 6.103785e-30 2.197362e-28


## Testing the effects of the habitat on community structure

vegan_otu <- function(physeq){
  OTU <- otu_table(physeq)
  if(taxa_are_rows(OTU)){
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

vegan_matrix_cryosphere <- vegan_otu(NOMIS_CRYOS_sub)
NOMISCRYO_bray <- vegdist(vegan_matrix_cryosphere, method="bray")
ado_habitat <- adonis2(NOMISCRYO_bray ~ Habitat, data=metadata_nmds_cryo, method="bray", permutations=999)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = NOMISCRYO_bray ~ Habitat, data = metadata_nmds_cryo, permutations = 999, method = "bray")
# Df SumOfSqs      R2      F Pr(>F)    
# Habitat    8   40.743 0.21161 13.387  0.001 ***
#   Residual 399  151.793 0.78839                  
# Total    407  192.535 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Pairwise adonis
library(pairwiseAdonis)
pairwise_habitat <- pairwise.adonis(NOMISCRYO_bray, metadata_nmds_cryo$Habitat, p.adjust.m="holm")
# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1          GFS vs Microbial mat  1  2.851253  7.210056 0.04500377   0.001      0.036   .
# 2               GFS vs Sediment  1  2.723501  7.005175 0.04578391   0.001      0.036   .
# 3                  GFS vs Water  1 14.387782 38.755060 0.15094176   0.001      0.036   .
# 4             GFS vs Permafrost  1  5.375110 13.129454 0.06626705   0.001      0.036   .
# 5                GFS vs Glacier  1  3.299011  8.613379 0.05396402   0.001      0.036   .
# 6                   GFS vs Snow  1  4.863377 11.998161 0.06416192   0.001      0.036   .
# 7            GFS vs Terrestrial  1  6.961109 17.954327 0.09024346   0.001      0.036   .
# 8             GFS vs Cryoconite  1  5.709281 14.299240 0.07474802   0.001      0.036   .
# 9     Microbial mat vs Sediment  1  2.126911  6.776995 0.28502318   0.001      0.036   .
# 10       Microbial mat vs Water  1  3.685020 11.115575 0.11102743   0.001      0.036   .
# 11  Microbial mat vs Permafrost  1  2.124891  4.894733 0.08038024   0.001      0.036   .
# 12     Microbial mat vs Glacier  1  2.782821  9.554946 0.30280343   0.001      0.036   .
# 13        Microbial mat vs Snow  1  2.068414  4.877341 0.09586471   0.001      0.036   .
# 14 Microbial mat vs Terrestrial  1  2.813035  7.802011 0.13046403   0.001      0.036   .
# 15  Microbial mat vs Cryoconite  1  2.335338  5.825020 0.10822141   0.001      0.036   .
# 16            Sediment vs Water  1  2.966834  9.443012 0.10326663   0.001      0.036   .
# 17       Sediment vs Permafrost  1  2.265084  5.396084 0.09919987   0.001      0.036   .
# 18          Sediment vs Glacier  1  2.707289 15.236753 0.50391499   0.001      0.036   .
# 19             Sediment vs Snow  1  2.250761  5.567717 0.12492713   0.001      0.036   .
# 20      Sediment vs Terrestrial  1  2.724959  8.171289 0.15367859   0.001      0.036   .
# 21       Sediment vs Cryoconite  1  2.425448  6.414965 0.13529411   0.001      0.036   .
# 22          Water vs Permafrost  1  6.860327 18.554571 0.13295566   0.001      0.036   .
# 23             Water vs Glacier  1  4.915089 15.934857 0.15480526   0.001      0.036   .
# 24                Water vs Snow  1  6.077502 16.892240 0.13208182   0.001      0.036   .
# 25         Water vs Terrestrial  1  8.956743 26.749233 0.18608261   0.001      0.036   .
# 26          Water vs Cryoconite  1  7.088575 20.190989 0.15159426   0.001      0.036   .
# 27        Permafrost vs Glacier  1  3.308816  8.256560 0.13262153   0.001      0.036   .
# 28           Permafrost vs Snow  1  2.418614  5.430222 0.06508699   0.001      0.036   .
# 29    Permafrost vs Terrestrial  1  4.269586 10.554020 0.11161895   0.001      0.036   .
# 30     Permafrost vs Cryoconite  1  3.288571  7.630748 0.08707843   0.001      0.036   .
# 31              Glacier vs Snow  1  3.251183  8.495850 0.16183851   0.001      0.036   .
# 32       Glacier vs Terrestrial  1  4.045746 12.581132 0.20103714   0.001      0.036   .
# 33        Glacier vs Cryoconite  1  3.504765  9.727343 0.17455243   0.001      0.036   .
# 34          Snow vs Terrestrial  1  3.933444  9.975410 0.11878965   0.001      0.036   .
# 35           Snow vs Cryoconite  1  3.044213  7.181098 0.09304218   0.001      0.036   .
# 36    Terrestrial vs Cryoconite  1  4.692128 12.332682 0.13961630   0.001      0.036   .
