###Computing the average ASV richness from rarefied data
library(dplyr)
library(tidyverse)
library(breakaway)
library(data.table)
library(phyloseq)
library(microViz)
library(viridis)
library(hrbrthemes)
library(paletteer)

NOMIS_FR <- readRDS("NOMIS_FR_20230901.RData")

uganda=c("Uganda")
prune_Uganda <- subset_samples(NOMIS_FR, !Site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))

## Calculate diversity metrics
diversity_nomis <- plot_richness(prune_Uganda, measures=c("Observed"), color="Site_c")
alphadt_nomis<- data.table(diversity_nomis$data)

shannon_gfs <- diversity(t(prune_Uganda_df), index = "shannon",base=2,71828)

## Plot ASV richness per region using violin plot and include jitter 
plot_div <- ggplot(alphadt_nomis,aes(x=alphadt_nomis$Site_c,y=value, color=Site_c)) + 
  geom_violin(scale="width") + 
  geom_jitter(aes(group=Site_c), position=position_jitterdodge()) + 
  stat_summary(fun.y="mean",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
               width=1, position=position_dodge(),show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, margin=margin(0.5, unit="cm"))) + theme_bw()

plot_div + scale_colour_manual(values=c("#2E2A2BFF","#CF4E9CFF","#8C57A2FF",
                               "#3EBCB6","#82581FFF","#2F509EFF",
                               "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E"))


## Compute the average observed richness per mountain range
## Filter asv richness
OR <- subset(alphadt_nomis, variable == "Observed")
hist(OR$value) ##normally distributed
ASVrichness_GFS <- OR %>% 
  group_by(Site_c) %>% 
  summarise(average=mean(value), std=sd(value))

# # A tibble: 10 × 3
# Site_c       average   std
# <chr>          <dbl> <dbl>
# 1 Alaska          959.  219.
# 2 Alps           1756.  497.
# 3 Caucasus       1820.  377.
# 4 Chile          1225.  559.
# 5 Ecuador        1549.  577.
# 6 Greenland      1163.  309.
# 7 kyrgyzstan   1383.  466.
# 8 Nepal          1873.  306.
# 9 New_Zealand    1619.  362.
# 10 Norway        1677.  667.
# 
average_richness <- alphadt_nomis %>% 
  summarise(average=mean(value), std=sd(value))
# average      std
# 1 1556.095 512.6902

median_richness <- alphadt_nomis %>% 
  summarise(median=median(value), x = quantile(value, c(0.25, 0.5, 0.75)))


Shannon <- subset(alphadt_nomis, variable == "Shannon")
hist(Shannon$value)
Shannon_GFS <- Shannon %>% 
  group_by(Site_c) %>% 
  summarise(average=mean(value), std=sd(value))
# 
# Site_c       average   std
# <chr>          <dbl> <dbl>
# 1 Alaska          5.42 0.692
# 2 Alps            6.13 0.468
# 3 Caucasus        6.01 0.494
# 4 Chile           5.39 0.906
# 5 Ecuador         5.94 0.544
# 6 Greenland       5.58 0.646
# 7 kyrgyzstan    5.85 0.458
# 8 Nepal           6.28 0.350
# 9 New_Zealand     5.97 0.522
# 10 Norway         6.00 0.668

median_shannon <- Shannon %>% 
  summarise(median=median(value), x = quantile(value, c(0.25, 0.5, 0.75)))


## Plot ASV richness per region using violin plot and include jitter 
plot_shannon <- ggplot(Shannon,aes(x=Site_c,y=value, color=Site_c)) + 
  geom_violin(scale="width") + 
  geom_jitter(aes(group=Site_c), position=position_jitterdodge()) + 
  stat_summary(fun.y="mean",geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), 
               width=1, position=position_dodge(),show.legend = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, margin=margin(0.5, unit="cm"))) + theme_bw()
plot_shannon + scale_colour_manual(values=c("#2E2A2BFF","#CF4E9CFF","#8C57A2FF",
                                        "#3EBCB6","#82581FFF","#2F509EFF",
                                        "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E"))

## Compute some statistics
## Evenness with Bulla index
bulla_estimate <- microbiome::evenness(prune_Uganda, index="all")
alphadt_nomis$even <-bulla_estimate$bulla
hist(alphadt_nomis$even)
## Compute the average evenness per mountain range
evenness_GFS <- alphadt_nomis %>% 
  group_by(Site_c) %>% 
  summarise(average=mean(even), std=sd(even))

evenness_GFS_median <- alphadt_nomis %>% 
  #summarise(average=mean(even), std=sd(even)) %>% 
  summarise(median=median(even), x = quantile(even, c(0.25, 0.5, 0.75)))
# 
# # A tibble: 10 × 3
# Site_c       average    std
# <chr>          <dbl>  <dbl>
#   1 Alaska         0.402 0.0852
# 2 Alps           0.414 0.0385
# 3 Caucasus       0.398 0.0454
# 4 Chile          0.376 0.0764
# 5 Ecuador        0.413 0.0472
# 6 Greenland      0.385 0.0743
# 7 kyrgyzstan   0.399 0.0442
# 8 Nepal          0.430 0.0383
# 9 New_Zealand    0.411 0.0498
# 10 Norway        0.413 0.0632

evenness_GFS_total <- alphadt_nomis %>% 
  summarise(average=mean(even), std=sd(even))

##lets do some stat
library(performance)
library(multcomp)
region_effect_anova_even <-aov(even ~ Site_c, data = alphadt_nomis)
summary(region_effect_anova_even)
aov_region_even <-anova(region_effect_anova_even)
check_model(region_effect_anova_even)
check_normality(region_effect_anova_even)
check_heteroscedasticity(region_effect_anova_even)

region_effect_kruskal <-kruskal.test(even ~ Site_c, data = alphadt_nomis)
summary(region_effect_kruskal)
# Kruskal-Wallis rank sum test
# 
# data:  even by Site_c
# Kruskal-Wallis chi-squared = 7.7787, df = 9, p-value = 0.5566

## Richness
region_effect_anova_richness <-aov(value ~ Site_c, data = OR)
summary(region_effect_anova_richness)
aov_region_richness <-anova(region_effect_anova_richness)
check_model(region_effect_anova_richness)
check_normality(region_effect_anova_richness)
check_heteroscedasticity(region_effect_anova_richness)

tukey_richness <- TukeyHSD(region_effect_anova_richness)
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = (value) ~ Site_c, data = OR)
# 
# $Site_c
# diff           lwr         upr     p adj
# Alps-Alaska               797.83333   340.8806463 1254.786020 0.0000046
# Caucasus-Alaska           860.91228   374.1368619 1347.687700 0.0000033
# Chile-Alaska              266.23333  -309.1217910  841.588458 0.8946102
# Ecuador-Alaska            590.77778    -3.4464395 1185.001995 0.0527045
# Greenland-Alaska          204.00000  -476.7693638  884.769364 0.9937561
# kyrgyzstan-Alaska       423.92157   -75.3266335  923.169771 0.1718806
# Nepal-Alaska              914.52083   408.0122181 1421.029449 0.0000019
# New_Zealand-Alaska        660.53333   179.1566998 1141.909967 0.0008338
# Norway-Alaska             718.66667   124.4424494 1312.890884 0.0058564
# Caucasus-Alps              63.07895  -362.2788448  488.436739 0.9999789
# Chile-Alps               -531.60000 -1056.0163060   -7.183694 0.0440823
# Ecuador-Alps             -207.05556  -752.1072625  337.996151 0.9676413
# Greenland-Alps           -593.83333 -1232.1325340   44.465867 0.0916576
# kyrgyzstan-Alps        -373.91176  -813.4885046   65.664975 0.1700355
# Nepal-Alps                116.68750  -331.1181724  564.493172 0.9978252
# New_Zealand-Alps         -137.30000  -556.4687086  281.868709 0.9881422
# Norway-Alps               -79.16667  -624.2183736  465.885040 0.9999823
# Chile-Caucasus           -594.67895 -1145.2759072  -44.081988 0.0232435
# Ecuador-Caucasus         -270.13450  -840.4203203  300.151314 0.8804443
# Greenland-Caucasus       -656.91228 -1316.8898711    3.065310 0.0521559
# kyrgyzstan-Caucasus    -436.99071  -907.4928486   33.511424 0.0928478
# Nepal-Caucasus             53.60855  -424.5906427  531.807748 0.9999981
# New_Zealand-Caucasus     -200.37895  -651.8730852  251.115190 0.9164983
# Norway-Caucasus          -142.24561  -712.5314314  428.040203 0.9984522
# Ecuador-Chile             324.54444  -322.9963837  972.085273 0.8401035
# Greenland-Chile           -62.23333  -790.0063959  665.539729 0.9999998
# kyrgyzstan-Chile        157.68824  -403.9659921  719.342463 0.9961562
# Nepal-Chile               648.28750    80.1698352 1216.405165 0.0123782
# New_Zealand-Chile         394.30000  -151.5297969  940.129797 0.3808358
# Norway-Chile              452.43333  -195.1074949 1099.974162 0.4301341
# Greenland-Ecuador        -386.77778 -1129.5580493  356.002494 0.8075338
# kyrgyzstan-Ecuador     -166.85621  -747.8246871  414.112269 0.9954421
# Nepal-Ecuador             323.74306  -263.4763092  910.962420 0.7506143
# New_Zealand-Ecuador        69.75556  -495.9290460  635.440157 0.9999957
# Norway-Ecuador            127.88889  -536.4739829  792.251761 0.9998090
# kyrgyzstan-Greenland    219.92157  -449.3084955  889.151633 0.9878675
# Nepal-Greenland           710.52083    35.8571477 1385.184519 0.0303083
# New_Zealand-Greenland     456.53333  -199.4724401 1112.539107 0.4360394
# Norway-Greenland          514.66667  -228.1136049 1257.446938 0.4425545
# Nepal-kyrgyzstan        490.59926    -0.2906714  981.489201 0.0502705
# New_Zealand-kyrgyzstan  236.61176  -228.3026526  701.526182 0.8275646
# Norway-kyrgyzstan       294.74510  -286.2233799  875.713576 0.8302011
# New_Zealand-Nepal        -253.98750  -726.6899702  218.714970 0.7774907
# Norway-Nepal             -195.85417  -783.0735314  391.365198 0.9865367
# Norway-New_Zealand         58.13333  -507.5512682  623.817935 0.9999991

## Shannon
region_effect_anova_shannon <-aov(sqrt(value) ~ Site_c, data = Shannon)
summary(region_effect_anova_shannon)
aov_region_shannon <-anova(region_effect_anova_shannon)
check_model(region_effect_anova_shannon)
check_normality(region_effect_anova_shannon)
check_heteroscedasticity(region_effect_anova_shannon)
##Ok we have to use non-parametric test

kruskal_shannon <-kruskal.test(value ~ Site_c, data = Shannon)
summary(kruskal_shannon)
library(FSA)
##Posthoc
DT_shannon = dunnTest(value ~ Site_c, data = Shannon,
              method="holm")     
#Results from Dunn Test
# Comparison           Z      P.unadj       P.adj
# 1               Alaska - Alps -3.31115485 9.291178e-04 0.040881183
# 2           Alaska - Caucasus -2.44953005 1.430428e-02 0.557866817
# 3             Alps - Caucasus  0.75388319 4.509194e-01 1.000000000
# 4              Alaska - Chile -0.27996795 7.795021e-01 1.000000000
# 5                Alps - Chile  2.57802836 9.936584e-03 0.407399956
# 6            Caucasus - Chile  1.87303980 6.106289e-02 1.000000000
# 7            Alaska - Ecuador -1.65122280 9.869309e-02 1.000000000
# 8              Alps - Ecuador  0.97576894 3.291790e-01 1.000000000
# 9          Caucasus - Ecuador  0.37029579 7.111621e-01 1.000000000
# 10            Chile - Ecuador -1.26650791 2.053313e-01 1.000000000
# 11         Alaska - Greenland -0.31764902 7.507512e-01 1.000000000
# 12           Alps - Greenland  2.03164187 4.218992e-02 1.000000000
# 13       Caucasus - Greenland  1.47902794 1.391328e-01 1.000000000
# 14          Chile - Greenland -0.07579935 9.395787e-01 0.939578733
# 15        Ecuador - Greenland  1.02984810 3.030813e-01 1.000000000
# 16      Alaska - kyrgyzstan -1.43410277 1.515429e-01 1.000000000
# 17        Alps - kyrgyzstan  1.81326217 6.979139e-02 1.000000000
# 18    Caucasus - kyrgyzstan  1.01253054 3.112845e-01 1.000000000
# 19       Chile - kyrgyzstan -0.98796058 3.231720e-01 1.000000000
# 20     Ecuador - kyrgyzstan  0.45651934 6.480166e-01 1.000000000
# 21   Greenland - kyrgyzstan -0.74672005 4.552326e-01 1.000000000
# 22             Alaska - Nepal -3.89168808 9.954918e-05 0.004479713
# 23               Alps - Nepal -1.02306081 3.062791e-01 1.000000000
# 24           Caucasus - Nepal -1.62861529 1.033945e-01 1.000000000
# 25              Chile - Nepal -3.18612262 1.441935e-03 0.062003192
# 26            Ecuador - Nepal -1.68587247 9.182037e-02 1.000000000
# 27          Greenland - Nepal -2.60118909 9.290123e-03 0.390185177
# 28       kyrgyzstan - Nepal -2.55698929 1.055825e-02 0.422329878
# 29       Alaska - New_Zealand -2.39051106 1.682494e-02 0.639347815
# 30         Alps - New_Zealand  0.86434157 3.874003e-01 1.000000000
# 31     Caucasus - New_Zealand  0.09221571 9.265267e-01 1.000000000
# 32        Chile - New_Zealand -1.81312045 6.981325e-02 1.000000000
# 33      Ecuador - New_Zealand -0.29970692 7.644007e-01 1.000000000
# 34    Greenland - New_Zealand -1.42451558 1.542973e-01 1.000000000
# 35 kyrgyzstan - New_Zealand -0.93514616 3.497130e-01 1.000000000
# 36        Nepal - New_Zealand  1.73563166 8.262898e-02 1.000000000
# 37            Alaska - Norway -1.81832556 6.901439e-02 1.000000000
# 38              Alps - Norway  0.79359081 4.274337e-01 1.000000000
# 39          Caucasus - Norway  0.19617871 8.444703e-01 1.000000000
# 40             Chile - Norway -1.41985191 1.556508e-01 1.000000000
# 41           Ecuador - Norway -0.14946125 8.811897e-01 1.000000000
# 42         Greenland - Norway -1.16353031 2.446144e-01 1.000000000
# 43      kyrgyzstan - Norway -0.62743482 5.303743e-01 1.000000000
# 44             Nepal - Norway  1.51677637 1.293232e-01 1.000000000
# 45       New_Zealand - Norway  0.12417358 9.011778e-01 1.000000000

##Results from Inext
#         Site         Diversity  Observed Estimator    s.e.       LCL       UCL
# 1       alps  Species richness 14130.000 22289.342 272.771 21771.999 22841.708
# 2       alps Shannon diversity  8760.813 11254.962  57.395 11142.470 11367.455
# 3       alps Simpson diversity  5835.264  6432.662  33.195  6367.600  6497.724
# 4         nz  Species richness 10802.000 15755.292 191.223 15394.460 16144.476
# 5         nz Shannon diversity  7256.587  9307.307  51.341  9206.682  9407.933
# 6         nz Simpson diversity  5140.499  5804.794  33.577  5738.983  5870.604
# 7   caucasus  Species richness 12069.000 19926.149 277.321 19401.142 20488.749
# 8   caucasus Shannon diversity  7841.318 10490.330  66.527 10359.939 10620.721
# 9   caucasus Simpson diversity  5363.305  6013.883  35.418  5944.465  6083.300
# 10        kh  Species richness 12332.000 20760.372 286.732 20216.865 21341.343
# 11        kh Shannon diversity  8616.254 12247.449  80.351 12089.964 12404.934
# 12        kh Simpson diversity  6056.276  7115.446  45.614  7026.044  7204.849
# 13 greenland  Species richness  4033.000  7871.895 217.448  7468.817  8322.261
# 14 greenland Shannon diversity  3307.025  5632.570  75.578  5484.440  5780.700
# 15 greenland Simpson diversity  2683.018  3632.849  43.958  3546.694  3719.004
# 16    norway  Species richness  7764.000 15547.691 324.568 14937.105 16210.251
# 17    norway Shannon diversity  5958.206  9767.697  96.040  9579.462  9955.931
# 18    norway Simpson diversity  4534.387  5760.727  56.236  5650.506  5870.948
# 19   ecuador  Species richness  7121.000 13113.421 256.910 12630.674 13638.465
# 20   ecuador Shannon diversity  5453.244  8649.617  79.402  8493.992  8805.241
# 21   ecuador Simpson diversity  4132.851  5220.711  46.054  5130.448  5310.975
# 22     nepal  Species richness  8515.000 14356.444 244.208 13897.082 14855.013
# 23     nepal Shannon diversity  5649.826  7670.129  53.929  7564.430  7775.829
# 24     nepal Simpson diversity  3945.986  4463.041  31.286  4401.721  4524.361
# 25    alaska  Species richness  7072.000 14988.881 342.959 14344.729 15690.086
# 26    alaska Shannon diversity  5028.972  8386.781  90.233  8209.928  8563.635
# 27    alaska Simpson diversity  3442.734  4224.420  47.574  4131.177  4317.662
# 28     chile  Species richness  7583.000 16084.586 341.272 15441.589 16780.193
# 29     chile Shannon diversity  6287.391 12201.955 142.930 11921.818 12482.092
# 30     chile Simpson diversity  5023.538  7947.849 106.146  7739.807  8155.892

##Results from INEXT - observed and estimated richness
results_inext = read.csv("results_INEXT_2022111.csv",sep=",",header=T)
results_inext <- as.data.frame(results_inext)

names(results_inext) <- c("Site","diversity","asv_obs","estimated","SE","LCL","UCL")

inext_res <- results_inext |> pivot_longer(cols = -c(Site, diversity, LCL, UCL), names_to = "Type") 
head(inext_res)

merged_pivot <- merge(results_inext, inext_res, by.x="Site",by.y="Site")

##Assign SD of Percent to 0
merged_pivot$SE[merged_pivot$Type == 'Observed'] = 0

merged_pivot <- as.data.frame(merged_pivot)
merged_pivot$SE <- as.numeric(merged_pivot$SE)
merged_pivot$value <- as.numeric(merged_pivot$value)

##Plot everything!
regional_div <- ggplot(merged_pivot, aes(x=Site, y = value, fill= Type)) + 
  geom_col(position="dodge") + 
  labs(y="Average")+ geom_errorbar(aes(ymin=(value - SE), ymax=(value + SE)),
                                   width=.2, colour="red", 
                                   position=position_dodge(.9))

regional_div + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

##With Facet
facet_regional <- ggplot(merged_pivot, aes( x = Site, y = value, fill=Type) ) + 
  geom_col(position="dodge") + labs(y="Number of ASVs")+ 
  geom_errorbar(aes(ymin=(value - SE), ymax=(value + SE)), width=.2, colour="red", 
                position=position_dodge(.9))+
  facet_wrap( ~ Site ) + 
  xlab("Mountain Ranges")  
p1 + theme_bw() + 
  theme( axis.text.x = element_text( angle = 90,  hjust = 1 ) )


## Effect of altitude on species richness
## On check comment la latitude et l'altitude sont lilées entre elles 
## de la latitude

lat_ele <- ggplot(data=alphadt_nomis,aes(x=lat_sp,y=ele_sp))+
   geom_point()+
   geom_smooth(method="gam", formula = y ~ s(x, bs = 'tp'))

 lat_ele + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# ## We test the effect of latitude on altitude via a spline
 library(mgcv)
 mod_lat <- gam(data=alphadt_nomis, formula = ele_sp ~ s(lat_sp, bs='tp'))
 summary(mod_lat)
 alphadt_nomis$ele_resids = mod_lat$residuals
 
 # Family: gaussian 
 # Link function: identity 
 # 
 # Formula:
 #   ele_sp ~ s(lat_sp, bs = "tp")
 # 
 # Parametric coefficients:
 #   Estimate Std. Error t value Pr(>|t|)    
 # (Intercept)  2453.68      37.14   66.07   <2e-16 ***
 #   ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 # 
 # Approximate significance of smooth terms:
 #   edf Ref.df     F p-value    
 # s(lat_sp) 7.435   8.13 158.7  <2e-16 ***
 #   ---
 #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 # R-sq.(adj) =  0.898   Deviance explained = 90.3%
 # GCV = 2.1507e+05  Scale est. = 2.0273e+05  n = 147
 # 


# ## We try now the inverse and investigate how altitude influences latitude
# mod_ele <- gam(data=alphadt_nomis, formula = lat_sp ~ s(ele_sp, bs='tp'))
# summary(mod_ele)
# alphadt_nomis$lat_resids = mod_ele$residuals

# ## On detrend en fonction de la latitude -- Robin advice
# ele_sp <- lm(data=alphadt_nomis, formula = value ~ ele_sp + lat_sp)
# summary(ele_sp)
# check_model(ele_sp)
# check_normality(ele_sp)
# check_heteroscedasticity(ele_sp)

# null_model_elevation <- glm(value ~ 1, data=alphadt_nomis)
# anova(ele_sp, null_model_elevation, test="F")
# 
# spe_ele_plot <- ggplot(alphadt_nomis, aes(x=ele_sp, y=value)) + geom_point() + geom_smooth(method='lm', color="black")
# spe_ele_plot + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
# 
# ele_sp$coefficients
# (Intercept)        ele_sp   abs(lat_sp) 
# 1921.50376520   -0.01085261   -7.87143852 

# ## We retrieve the residuals of the models and investigate the correlations between altitude or latitude and species richness
pred <- predict(lm(data=alphadt_nomis, formula = value ~ ele_resids), 
                 se.fit = TRUE, interval = "confidence")
limits <- as.data.frame(pred$fit)
spe_ele <- ggplot(alphadt_nomis, aes(x=ele_resids, y=value)) + geom_point(size=3, alpha=0.4,color="#3F459BFF") + 
geom_smooth(method='gam', se=T, formula = y ~ s(x, bs='tp'), fill="blue", color="black", alpha=0.2, span=0.3)  +
geom_line(aes(x = ele_resids, y = limits$lwr), 
          linetype = 2) +
  geom_line(aes(x = ele_resids, y = limits$upr), 
            linetype = 2) 
spe_ele + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
 
model_elevation<-gam(data=alphadt_nomis, formula = value ~ s(ele_resids, bs='tp'))
summary(model_elevation)
model_elevation$coefficients

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   value ~ s(ele_resids, bs = "tp")
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1556.10      40.64   38.29   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(ele_resids)   1      1 13.09 0.00041 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0765   Deviance explained = 8.28%
# GCV = 2.461e+05  Scale est. = 2.4275e+05  n = 147

catchment_area <-read.csv("catchment_final.csv",sep=",",header=T)
merge_data <- merge(catchment_area, alphadt_nomis, by.x="Glacier", by.y="code_gl")

##There is no relationship between the %of glacier coverage and the altitude
cov_ele <- ggplot(data=merge_data,aes(x=(gl_cov),y=ele_sp))+
  geom_point()+
  geom_smooth(method="gam", formula = y ~ s(x, bs = 'tp'))

cov_ele + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# ## We test the effect % of glacier coverage on species diversity

pred_cov <- predict(lm(data=merge_data, formula = value ~ gl_cov), 
                se.fit = TRUE, interval = "confidence")
limits_cov<- as.data.frame(pred_cov$fit)

spe_cov <- ggplot(merge_data, aes(x=(gl_cov), y=value)) + geom_point(size=3, alpha=0.4,color="#3F459BFF") + 
  geom_smooth(method='gam', se=T, formula = y ~ s(x, bs='tp'), fill="blue", color="black", alpha=0.2, span=0.3)  +
  geom_line(aes(x = gl_cov, y = limits_cov$lwr), 
            linetype = 2) +
  geom_line(aes(x = gl_cov, y = limits_cov$upr), 
            linetype = 2) 
spe_cov + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

model_cov<-gam(data=merge_data, formula = value ~ s(gl_cov, bs='tp'))
summary(model_cov)
model_cov$coefficients
# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   value ~ s(gl_cov, bs = "tp")
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1556.10      40.73    38.2   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(gl_cov)   1      1 12.35 0.00059 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.0721   Deviance explained = 7.85%
# GCV = 2.4726e+05  Scale est. = 2.4389e+05  n = 147


