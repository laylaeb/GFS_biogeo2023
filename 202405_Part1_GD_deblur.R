library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(knitr)
library(dplyr)
library(data.table)
library(plyr)
library(DescTools)
library(devtools)
library(phyloseqCompanion)
library(iNEXT)

## Load ASV and taxonomy tables
asv<-read.csv(file="20240221_NOMIS_rarefied_deblur_table.csv.gz",sep=",",header=TRUE,row.names=1)
tax<-read.csv(file="20240221_NOMIS_rarefied_deblur_taxonomy.csv.gz",sep=",",header=TRUE,row.names=1)

## Load metadata file
metadata_NOMIS="/Users/ezzat/Leilas_Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_August2023/DEBLUR/revision_2/202402_NOMIS_metadata_GFS.tsv"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

## create OTU table 
OTU_NOMIS <- otu_table(asv, taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))

## Create phyloseq object
merged_NOMIS_DEBLUR <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

## Prune Uganda samples
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR, !site_c %in% "Uganda")
metadata_nomis_inext <- sample.data.frame(prune_Uganda)

sites <- c("Greenland", "Norway", "Alps", "New_Zealand", "Ecuador", "Caucasus", "Kirghizistan", "Chile", "Alaska", "Nepal")

sample_lists <- list()
asv_lists <- list()
vec_lists <- list()

# Loop through each site and perform the subset and merge operations
for (site in sites) {
   subset_result <- subset_samples(prune_Uganda, site_c == site)
   merged_result <- merge_samples(subset_result,"sample")

   # ## create asv table
   asv_table_site <- otu_table(merged_result, taxa_are_rows=T)
   asv_table_site_t <- t(asv_table_site)
   asv_lists[[site]] <- asv_table_site_t
   #
   # ## calculate rowSums and
   asv_table_df <- asv_table_site_t
   sumrow_site <- unname(rowSums(asv_table_df>0))
   sort_site<- sort(sumrow_site, decreasing=T)
   vec_site <- sort_site[sort_site >0]
   vec_lists[[site]] <- vec_site
}

list_exped_all <- list(alps=c(ncol(asv_lists$Alps),vec_lists$Alps),nz=c(ncol(asv_lists$New_Zealand),vec_lists$New_Zealand),
                       caucasus=c(ncol(asv_lists$Caucasus),vec_lists$Caucasus),
                       kh=c(ncol(asv_lists$Kirghizistan),vec_lists$Kirghizistan),
                       greenland=c(ncol(asv_lists$Greenland),vec_lists$Greenland),
                       norway=c(ncol(asv_lists$Norway),vec_lists$Norway), 
                       ecuador=c(ncol(asv_lists$Ecuador),vec_lists$Ecuador),
                       nepal=c(ncol(asv_lists$Nepal),vec_lists$Nepal),
                       alaska=c(ncol(asv_lists$Alaska),vec_lists$Alaska),
                       chile=c(ncol(asv_lists$Chile),vec_lists$Chile))
                       

out_all_exped <- iNEXT(list_exped_all, q=0, datatype="incidence_freq", se=T, conf=0.95, nboot=99)

# Compare 10 assemblages with Hill number order q = 0.
# $class: iNEXT
# 
# $DataInfo: basic data information
# Assemblage  T     U S.obs     SC   Q1   Q2   Q3   Q4   Q5  Q6  Q7  Q8  Q9 Q10
# 1        alps 26 74090 14345 0.9553 3524 2789 1662 1046  799 632 546 423 382 344
# 2          nz 20 62266 12839 0.9749 1779 2364 2073 1607 1112 798 613 508 379 323
# 3    caucasus 19 69809 13679 0.9702 2334 2541 1976 1349 1006 754 603 521 391 362
# 4          kh 18 48162 11138 0.9433 2966 2165 1469  952  660 490 405 327 290 283
# 5   greenland  7 18457  8594 0.7752 4574 1405  997  665  451 347 155   0   0   0
# 6      norway  9 22731  7146 0.9417 1675 1766 1332  792  522 356 254 205 244   0
# 7     ecuador 10 20338  7847 0.8715 3000 1995 1025  663  390 260 185 172 130  27
# 8       nepal 17 49491 13174 0.9391 3313 2642 2105 1475  964 671 445 363 298 218
# 9      alaska 15 35312 10110 0.9131 3326 1947 1321  900  641 475 350 272 217 181
# 10      chile 10 21176  9865 0.8366 4046 3077 1392  623  368 178  97  46  27  11
# 
# $iNextEst: diversity estimates with rarefied and extrapolated samples.
# $size_based (LCL and UCL are obtained for fixed size.)
# 
# Assemblage  t        Method Order.q        qD    qD.LCL    qD.UCL        SC    SC.LCL    SC.UCL
# 1         alps  1   Rarefaction       0  2849.615  2832.341  2866.889 0.3714909 0.3690710 0.3739108
# 10        alps 13   Rarefaction       0 11658.568 11593.878 11723.259 0.8928471 0.8915116 0.8941825
# 20        alps 26      Observed       0 14345.000 14259.267 14430.733 0.9552684 0.9539001 0.9566367
# 30        alps 38 Extrapolation       0 15460.967 15349.536 15572.397 0.9785872 0.9772830 0.9798913
# 40        alps 52 Extrapolation       0 16051.854 15908.561 16195.148 0.9909341 0.9899940 0.9918742
# 41          nz  1   Rarefaction       0  3113.300  3092.763  3133.837 0.3628507 0.3600469 0.3656546
# 50          nz 10   Rarefaction       0 11078.189 11027.705 11128.673 0.9014820 0.9000220 0.9029420
# 60          nz 20      Observed       0 12839.000 12777.272 12900.728 0.9749351 0.9735964 0.9762737
# 70          nz 30 Extrapolation       0 13303.195 13221.088 13385.303 0.9932316 0.9921647 0.9942985
# 80          nz 40 Extrapolation       0 13428.544 13327.581 13529.506 0.9981723 0.9975473 0.9987974
# 81    caucasus  1   Rarefaction       0  3674.158  3653.362  3694.954 0.4326679 0.4300568 0.4352789
# 90    caucasus 10   Rarefaction       0 11809.594 11754.847 11864.341 0.9111866 0.9100578 0.9123155
# 99    caucasus 19      Observed       0 13679.000 13612.260 13745.740 0.9701738 0.9690413 0.9713064
# 109   caucasus 29 Extrapolation       0 14370.351 14285.830 14454.871 0.9904792 0.9894736 0.9914848
# 118   caucasus 38 Extrapolation       0 14578.521 14476.199 14680.842 0.9965932 0.9959298 0.9972566
# 119         kh  1   Rarefaction       0  2675.667  2657.324  2694.009 0.4115009 0.4080310 0.4149708
# 127         kh  9   Rarefaction       0  8942.713  8889.529  8995.897 0.8680322 0.8662252 0.8698392
# 136         kh 18      Observed       0 11138.000 11069.228 11206.772 0.9432865 0.9414953 0.9450776
# 145         kh 27 Extrapolation       0 12142.667 12052.606 12232.729 0.9729811 0.9712557 0.9747065
# 154         kh 36 Extrapolation       0 12621.301 12506.224 12736.378 0.9871279 0.9858468 0.9884091
# 155  greenland  1   Rarefaction       0  2636.714  2606.333  2667.095 0.3856893 0.3803644 0.3910142
# 158  greenland  4   Rarefaction       0  6404.514  6332.876  6476.152 0.6906287 0.6853779 0.6958795
# 161  greenland  7      Observed       0  8594.000  8489.770  8698.230 0.7751983 0.7695096 0.7808869
# 165  greenland 11 Extrapolation       0 10654.608 10506.987 10802.230 0.8477848 0.8417159 0.8538536
# 168  greenland 14 Extrapolation       0 11750.283 11567.812 11932.754 0.8863808 0.8807471 0.8920144
# 169     norway  1   Rarefaction       0  2525.667  2502.595  2548.739 0.4501782 0.4458290 0.4545273
# 173     norway  5   Rarefaction       0  6037.508  5989.815  6085.201 0.8467197 0.8432725 0.8501669
# 177     norway  9      Observed       0  7146.000  7087.042  7204.958 0.9416833 0.9385746 0.9447921
# 182     norway 14 Extrapolation       0  7632.885  7555.457  7710.312 0.9818960 0.9795823 0.9842096
# 186     norway 18 Extrapolation       0  7766.099  7675.264  7856.934 0.9928983 0.9914098 0.9943868
# 187    ecuador  1   Rarefaction       0  2033.800  2010.194  2057.406 0.3435714 0.3387148 0.3484281
# 191    ecuador  5   Rarefaction       0  5800.917  5745.211  5856.623 0.7331303 0.7286513 0.7376093
# 196    ecuador 10      Observed       0  7847.000  7771.647  7922.353 0.8714846 0.8667721 0.8761971
# 201    ecuador 15 Extrapolation       0  8857.960  8758.302  8957.619 0.9354842 0.9307943 0.9401741
# 206    ecuador 20 Extrapolation       0  9365.471  9235.984  9494.958 0.9676125 0.9639889 0.9712362
# 207      nepal  1   Rarefaction       0  2911.235  2891.153  2931.318 0.3363086 0.3336077 0.3390095
# 215      nepal  9   Rarefaction       0 10843.825 10786.946 10900.704 0.8542432 0.8522143 0.8562721
# 223      nepal 17      Observed       0 13174.000 13099.927 13248.073 0.9391266 0.9370409 0.9412123
# 232      nepal 26 Extrapolation       0 14297.747 14192.627 14402.867 0.9741167 0.9722767 0.9759567
# 240      nepal 34 Extrapolation       0 14740.330 14607.219 14873.441 0.9878974 0.9866073 0.9891874
# 241     alaska  1   Rarefaction       0  2354.133  2335.859  2372.408 0.3712498 0.3671945 0.3753050
# 248     alaska  8   Rarefaction       0  8038.573  7975.991  8101.155 0.8328249 0.8305557 0.8350942
# 255     alaska 15      Observed       0 10110.000 10029.298 10190.702 0.9130799 0.9107569 0.9154029
# 263     alaska 23 Extrapolation       0 11366.869 11265.014 11468.724 0.9542825 0.9519046 0.9566605
# 270     alaska 30 Extrapolation       0 11966.604 11842.254 12090.954 0.9739430 0.9720315 0.9758546
# 271      chile  1   Rarefaction       0  2117.600  2093.860  2141.340 0.2322650 0.2279251 0.2366049
# 275      chile  5   Rarefaction       0  7025.929  6959.369  7092.488 0.6406129 0.6356874 0.6455384
# 280      chile 10      Observed       0  9865.000  9768.232  9961.768 0.8365567 0.8316015 0.8415119
# 285      chile 15 Extrapolation       0 11162.434 11036.033 11288.835 0.9251326 0.9207719 0.9294933
# 290      chile 20 Extrapolation       0 11756.741 11602.269 11911.212 0.9657060 0.9625420 0.9688699
# 
# NOTE: The above output only shows five estimates for each assemblage; call iNEXT.object$iNextEst$size_based to view complete output.
# 
# $coverage_based (LCL and UCL are obtained for fixed coverage; interval length is wider due to varying size in bootstraps.)
# 
# Assemblage        SC  t        Method Order.q        qD    qD.LCL    qD.UCL
# 1         alps 0.3714909  1   Rarefaction       0  2849.615  2832.341  2866.889
# 10        alps 0.8928470 13   Rarefaction       0 11658.565 11577.694 11739.436
# 20        alps 0.9552684 26      Observed       0 14345.000 14221.209 14468.791
# 30        alps 0.9785872 38 Extrapolation       0 15460.967 15302.470 15619.463
# 40        alps 0.9909341 52 Extrapolation       0 16051.854 15869.558 16234.150
# 41          nz 0.3628507  1   Rarefaction       0  3113.300  3092.763  3133.837
# 50          nz 0.9014818 10   Rarefaction       0 11078.184 11022.338 11134.030
# 60          nz 0.9749351 20      Observed       0 12839.000 12749.259 12928.741
# 70          nz 0.9932316 30 Extrapolation       0 13303.195 13190.044 13416.346
# 80          nz 0.9981723 40 Extrapolation       0 13428.544 13308.097 13548.990
# 81    caucasus 0.4326679  1   Rarefaction       0  3674.158  3653.362  3694.954
# 90    caucasus 0.9111859 10   Rarefaction       0 11809.571 11743.066 11876.075
# 99    caucasus 0.9701738 19      Observed       0 13679.000 13511.961 13846.039
# 109   caucasus 0.9904792 29 Extrapolation       0 14370.351 14254.826 14485.876
# 118   caucasus 0.9965932 38 Extrapolation       0 14578.521 14453.196 14703.845
# 119         kh 0.4115009  1   Rarefaction       0  2675.667  2657.324  2694.009
# 127         kh 0.8680317  9   Rarefaction       0  8942.701  8873.058  9012.344
# 136         kh 0.9432865 18      Observed       0 11138.000 11020.887 11255.113
# 145         kh 0.9729811 27 Extrapolation       0 12142.667 12013.464 12271.870
# 154         kh 0.9871279 36 Extrapolation       0 12621.301 12472.901 12769.700
# 155  greenland 0.3856893  1   Rarefaction       0  2636.714  2606.333  2667.095
# 158  greenland 0.6906259  4   Rarefaction       0  6404.451  6279.154  6529.748
# 161  greenland 0.7751983  7      Observed       0  8594.000  8308.885  8879.115
# 165  greenland 0.8477848 11 Extrapolation       0 10654.608 10436.656 10872.561
# 168  greenland 0.8863808 14 Extrapolation       0 11750.283 11499.726 12000.840
# 169     norway 0.4501782  1   Rarefaction       0  2525.667  2502.595  2548.739
# 173     norway 0.8467184  5   Rarefaction       0  6037.492  5979.630  6095.354
# 177     norway 0.9416833  9      Observed       0  7146.000  7065.139  7226.861
# 182     norway 0.9818960 14 Extrapolation       0  7632.885  7531.934  7733.835
# 186     norway 0.9928983 18 Extrapolation       0  7766.099  7658.696  7873.501
# 187    ecuador 0.3435714  1   Rarefaction       0  2033.800  2010.194  2057.406
# 191    ecuador 0.7331292  5   Rarefaction       0  5800.901  5733.521  5868.280
# 196    ecuador 0.8714846 10      Observed       0  7847.000  7688.318  8005.682
# 201    ecuador 0.9354842 15 Extrapolation       0  8857.960  8712.652  9003.269
# 206    ecuador 0.9676125 20 Extrapolation       0  9365.471  9193.868  9537.074
# 207      nepal 0.3363086  1   Rarefaction       0  2911.235  2891.153  2931.318
# 215      nepal 0.8542423  9   Rarefaction       0 10843.803 10772.265 10915.340
# 223      nepal 0.9391266 17      Observed       0 13174.000 13071.832 13276.168
# 232      nepal 0.9741167 26 Extrapolation       0 14297.747 14148.760 14446.734
# 240      nepal 0.9878974 34 Extrapolation       0 14740.330 14573.528 14907.132
# 241     alaska 0.3712498  1   Rarefaction       0  2354.133  2335.859  2372.408
# 248     alaska 0.8328249  8   Rarefaction       0  8038.573  7953.168  8123.977
# 255     alaska 0.9130799 15      Observed       0 10110.000 10014.984 10205.016
# 263     alaska 0.9542825 23 Extrapolation       0 11366.869 11228.874 11504.864
# 270     alaska 0.9739430 30 Extrapolation       0 11966.604 11808.700 12124.508
# 271      chile 0.2322650  1   Rarefaction       0  2117.600  2093.860  2141.340
# 275      chile 0.6406118  5   Rarefaction       0  7025.912  6941.805  7110.018
# 280      chile 0.8365567 10      Observed       0  9865.000  9664.166 10065.834
# 285      chile 0.9251326 15 Extrapolation       0 11162.434 10991.470 11333.397
# 290      chile 0.9657060 20 Extrapolation       0 11756.741 11564.793 11948.688
# 
# NOTE: The above output only shows five estimates for each assemblage; call iNEXT.object$iNextEst$coverage_based to view complete output.
# 
# $AsyEst: asymptotic diversity estimates along with related statistics.
# Assemblage         Diversity  Observed Estimator      s.e.       LCL       UCL
# 1      alaska  Species richness 10110.000 12761.461 121.63244 12523.066 12999.857
# 2      alaska Shannon diversity  7323.372  8680.045  37.74688  8606.063  8754.028
# 3      alaska Simpson diversity  5697.784  6341.104  25.24304  6291.629  6390.580
# 4        alps  Species richness 14345.000 16485.720  99.95291 16289.816 16681.625
# 5        alps Shannon diversity  9502.136 10546.720  25.52659 10496.688 10596.751
# 6        alps Simpson diversity  7202.104  7670.754  22.02076  7627.594  7713.914
# 7    caucasus  Species richness 13679.000 14694.514  68.42342 14560.407 14828.622
# 8    caucasus Shannon diversity  9969.948 10891.656  26.19916 10840.307 10943.005
# 9    caucasus Simpson diversity  7943.654  8491.867  25.84167  8441.218  8542.516
# 10      chile  Species richness  9865.000 12259.070 102.76892 12057.646 12460.493
# 11      chile Shannon diversity  8213.036 10925.080  56.29732 10814.740 11035.421
# 12      chile Simpson diversity  6852.220  9117.169  52.25107  9014.759  9219.579
# 13    ecuador  Species richness  7847.000  9877.075 101.02540  9679.069 10075.081
# 14    ecuador Shannon diversity  6149.224  7725.678  49.20424  7629.240  7822.117
# 15    ecuador Simpson diversity  4970.012  5919.584  40.60090  5840.008  5999.160
# 16  greenland  Species richness  8594.000 14975.742 162.55280 14657.144 15294.339
# 17  greenland Shannon diversity  6813.804  9899.642  69.20969  9763.993 10035.290
# 18  greenland Simpson diversity  5569.175  6836.369  49.93229  6738.504  6934.234
# 19         kh  Species richness 11138.000 13056.805  95.48701 12869.654 13243.956
# 20         kh Shannon diversity  7788.571  8812.076  30.26740  8752.753  8871.399
# 21         kh Simpson diversity  6023.627  6502.213  24.64965  6453.901  6550.525
# 22      nepal  Species richness 13174.000 15129.020  97.22625 14938.460 15319.580
# 23      nepal Shannon diversity  9844.280 11335.793  37.96731 11261.379 11410.208
# 24      nepal Simpson diversity  7756.070  8656.439  33.40887  8590.959  8721.919
# 25     norway  Species richness  7146.000  7852.084  55.56177  7743.185  7960.983
# 26     norway Shannon diversity  5824.175  6675.263  26.74034  6622.853  6727.673
# 27     norway Simpson diversity  4939.991  5610.371  27.36369  5556.739  5664.002
# 28         nz  Species richness 12839.000 13474.913  61.26178 13354.843 13594.984
# 29         nz Shannon diversity  9781.128 10723.386  28.32102 10667.878 10778.894
# 30         nz Simpson diversity  7887.599  8580.112  28.35501  8524.537  8635.687

df <- fortify(out_all_exped, type =1)

df.point <- df[which(df$Method=="Observed"),]
df.line <- df[which(df$Method!="Observed"),]
df.line$Method <- factor(df.line$Method, 
                         c("Rarefaction", "Extrapolation"),
                       )

df.asympote <- data.frame(y = c(24,10),
                          Asymptote = c("alps","nz","caucasus","kh","greenland","norway","ecuador","nepal","alaska","chile"))


ggplot(df, aes(x=x, y=y, colour=Assemblage)) + 
  #geom_point(aes(shape=Assemblage), size=5, data=df.point) +
  geom_line(aes(linetype= Method), lwd=1.5, data=df.line) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=Assemblage, colour=NULL), alpha=0.2) +
  labs(x="Number of GFS", y="Species diversity") +
scale_fill_manual(values=c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
                              "#3EBCB6","#82581FFF","#2F509EFF",
                              "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
)+
  scale_color_manual(values=c("#2E2A2BFF", "#CF4E9CFF","#8C57A2FF",
                             "#3EBCB6","#82581FFF","#2F509EFF",
                             "#E5614CFF","#97A1A7FF","#DC9445FF","#bee183")
  )+
  scale_linetype_discrete(name ="Method")+
theme_bw() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


inext_freq_results <- out_all_exped$AsyEst  
inext_freq_results$prop <- inext_freq_results$Observed/inext_freq_results$Estimator
inext_freq_results<-inext_freq_results[inext_freq_results$Diversity == 'Species richness',]
median_GD_freq <- inext_freq_results %>% 
  summarise(med = median(prop), 
            lower_quartile = quantile(prop, 0.25),
            median = quantile(prop, 0.5),
            upper_quartile = quantile(prop, 0.75))
# > median_GD_freq
# median         x
# 1 0.7361138 0.7260968
# 2 0.7361138 0.7361138
# 3 0.7361138 0.8441630

# med lower_quartile    median upper_quartile
# 1 0.8615944      0.7970271 0.8615944      0.9002519
# 
