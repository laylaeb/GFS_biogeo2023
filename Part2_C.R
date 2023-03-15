library("reshape2")
library("speedyseq")
library("phyloseq")
library("phyloseqCompanion")

###Core microbiome across the different regions
core_data <- as.data.frame(sample_data(prune_Uganda))
asv_count_level <- otu_table(prune_Uganda, taxa_are_rows=T)
asv_count_levelf <- asv_count_level[rowSums(asv_count_level[])>0,]

asv_core_ab  = transform_sample_counts(prune_Uganda, function(x) x / sum(x))
asv_core <- otu_table(asv_core_ab, taxa_are_rows=T)
asv_coref <- asv_core[rowSums(asv_core[])>0,]

tax_core <-as.matrix(tax_table(prune_Uganda))
write.csv(asv_core, "asv_core.csv")

metadata_nomis <- sample.data.frame(prune_Uganda)

# Calculating the binomial for the proba of occurrence
prob_occ <- function(asv_table, metadata, asv, ab_thr){
  subset = as.data.frame(t(as.matrix(asv_table[rownames(asv_table) == asv,])))
  
  subset[subset < ab_thr] = 0
  subset[subset > 0] = 1
  subset$Site_c = vapply(1:nrow(subset), function(x){as.character(metadata$Site_c[metadata$code_gl == rownames(subset)[x]])}, FUN.VALUE = character(1))
  colnames(subset)[colnames(subset) == asv] = 'asv'
  fit = glm(asv ~ Site_c, data=subset, family = binomial())
  # fit = gam(asv ~ s(latitude, longitude, type='sos'), data=subset, family=binomial())
  # pred_df = metadata[,c('latitude','longitude')]
  pred_df = expand.grid(Site_c=unique(metadata$Site_c))
  pred_df$pred = predict(fit, newdata=pred_df, type='response')
    
  pred_out = data.frame(ASV = asv, Ab_thr = ab_thr)
  pred_out$average_pred = mean(pred_df$pred)
  for (site_c in unique(metadata$Site_c)){
    pred_out[site_c] = pred_df$pred[pred_df$Site_c == site_c]}
  return(pred_out)}

seq_log <- function(v1, v2, n) {exp(seq(from = log(v1), to = log(v2), length.out = n))}

library(foreach)
library(doMC)
registerDoMC(4)
proba_df = foreach(abundance=seq_log(1/5000,1,30), .combine='rbind') %:% foreach(asv=rownames(asv_coref), .combine='rbind') %dopar% {prob_occ(asv_coref, core_data, asv, abundance)}

ggplot(data=proba_df, aes(x=(Ab_thr), y=(average_pred), group=ASV))+
  geom_line() + scale_x_log10() 

core_size <- expand.grid(Abundance=seq_log(1/5000,1,30),Prevalence+scale_x_log10()=seq_log(1/5000,1,30))
core_size$N = vapply(1:nrow(core_size), function(x){length(unique(proba_df$ASV[(proba_df$Ab_thr > core_size$Abundance[x]) & (proba_df$average_pred > core_size$Prevalence[x])]))}, FUN.VALUE = numeric(1))

ggplot(core_size, aes(x=Abundance,y=Prevalence)) + geom_tile(aes(fill=N)) +
  geom_line(data=data.frame(x=c(0,1),y=c(0,1)), aes(x=x,y=y), color='darkgrey', linetype = "dashed") + geom_point(aes(x=0.0008685272, y=0.2), color='red',size=5) +
  scale_x_log10() + scale_y_log10() + scale_fill_gradient(name = 'Core size', trans = "log", breaks = c(1,10,100,1000),na.value = 'transparent') + theme_linedraw()

ggsave('core_microbiome_biogeo.pdf', width=6, height=5)

##We can choose a threshold -- we would want the ASVs to be present in at least 5 out of 10 regions
write.csv(proba_df,"proba_df_20221026.csv")
proba_df_core <- read.csv("proba_df_20221026.csv", sep=",")
##Look at the abundance threshold of 0,1%
proba_df_thre <- 
  filter(proba_df_core, Ab_thr=="0.001165022")

##Then we want to keep taxa that are present at least in 5 regions out of 10! 
##In the Alps, there are 26 glaciers, so we mainly divided 1/26 to get the threshold 
binary_transform <- sapply(proba_df_thre[,4:13], function(x) ifelse(x > 0.03846154, TRUE, FALSE),
                          USE.NAMES = F)
binary_transform_merge <- merge(x=binary_transform, y=proba_df_thre, by="row.names")

binary_sum <- binary_transform_merge %>% 
  rowwise()%>% mutate(sum=sum(c_across(New_Zealand.x:Norway.x)))%>%filter(sum>5)

binary_tax_merge <- merge(x=binary_sum, y=tax_core, by.x="ASV", by.y=0)
write.csv(binary_tax_merge,"binary_tax_merge.csv")
#save(binary_tax_merge, file="binary_tax_merge.RData")
load("binary_tax_merge.RData")

##Relative abundance of the core microbiome
merge_core_abondance <- merge(asv_coref, binary_tax_merge, by.x="row.names", by.y="ASV" )
merge_core_ab <- read.csv("merge_NOMIS_core_asv.csv",sep=",",row.names=1)
merge_core_ab <- otu_table(merge_core_ab, taxa_are_rows=T)
tax_NOMIS <- tax_table(tax_core)

## Merge everything (mapping file, tree and OTUs) 
merged_NOMIS_core_ab<- merge_phyloseq(merge_core_ab, tax_NOMIS, metadata_nomis)
merged_NOMIS_core_ab_otu <- (as.matrix(otu_table(merged_NOMIS_core_ab, taxa_are_rows=T)))

melt_asv <- melt(merged_NOMIS_core_ab_otu)
merge_asv_data <- merge(as.data.frame(melt_asv),as.matrix(metadata_nomis), by.x="Var2",by.y="code_gl")

##Here we would need to divide by the number of glaciers per mountain ranges. 
sum_mr <- merge_asv_data %>% group_by(Site_c)%>% summarize(summrr=sum(value), n=n_distinct(Var2))%>% summarize(ar_mr=summrr/n, Site_c)
median(sum_mr$ar_mr)
quantile(sum_mr$ar_mr, prob=c(.25,.5,.75), type=1)
# # A tibble: 10 × 2
# ar_mr Site_c      
# <dbl> <chr>       
# 1 0.333 Alaska      
# 2 0.352 Alps        
# 3 0.346 Caucasus    
# 4 0.218 Chile       
# 5 0.138 Ecuador     
# 6 0.300 Greenland   
# 7 0.203 kyrgyzstan
# 8 0.295 Nepal       
# 9 0.180 New_Zealand 
# 10 0.278 Norway      


# > mean(sum_mr$ar_mr)
# [1] 0.2643127


#### CoreBarplot // what are the most abundant phlya 
### Families and Genera

asv_ab <- otu_table(merge_core_ab, taxa_are_rows=T)
merged_NOMIS_core_asv <- merge_phyloseq(asv_ab, tax_core, core_data)
transf_tax  = transform_sample_counts(merged_NOMIS_core_asv, function(x) x / sum(x))
sample_merge_region = merge_samples(transf_tax, "Site_c")
trans_ra = transform_sample_counts(sample_merge_region, function(x) x / sum(x))

##Turn all ASVs into Phylum counts
data <- psmelt(trans_ra) # create dataframe from phyloseq object
data$Phylum<- as.character(data$Phylum) #convert to character

##We kept the 15 most abundant genera
sumtot <- data %>% group_by(Phylum) %>% summarize(sum=sum(Abundance))%>%
  filter(!(Phylum %in% c("","g__uncultured"))) %>%
  slice_max(n=15, order_by=sum)

data$Phylum[!(data$Phylum %in% sumtot$Phylum)] <- "Other"

##Turn all ASVs into Family counts
data <- psmelt(trans_ra) # create dataframe from phyloseq object
data$Family<- as.character(data$Family) #convert to character

##We kept the 15 most abundant genera
sumtot <- data %>% group_by(Family) %>% summarize(sum=sum(Abundance))%>%
  filter(!(Family %in% c("","g__uncultured"))) %>%
  slice_max(n=15, order_by=sum)

data$Family[!(data$Family %in% sumtot$Family)] <- "Other"

##Generating n colors via rbrewer
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


barplot_biogeo <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=Family))
barplot_biogeo + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E",
                               "#E6AB02","#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                               "#FFFF99")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))


##What are the most abundant phlya?
### First we tax_glom the dataset
taxglom_biogeo_phylum <- tax_glom(physeq=merged_NOMIS_core_asv, taxrank=rank_names(merged_NOMIS_core_asv)[2], NArm=F)
tax_table_core_phylum <- tax_table(taxglom_biogeo_phylum)
## Take the relative abundance and group by mountain ranges
transf_ab_phy = transform_sample_counts(taxglom_biogeo_phylum, function(x) x / sum(x))
sample_merge_phy = merge_samples(transf_ab_phy, "Site_c")
trans_ra_phy = transform_sample_counts(sample_merge_phy, function(x) x / sum(x))

asv_core_phy <- otu_table(trans_ra_phy, taxa_are_rows=T)
melt_core_phy <- psmelt(asv_core_phy)
merge_taxo_core <- merge(melt_core_phy, tax_table_core_phylum, by.x="OTU", by.y="row.names")

##Now we would like to know what are the most abundant phyla across the dataset.
sumtot_phyla <- merge_taxo_core  %>% group_by(Phylum) %>% summarize(sum=sum(Abundance/10))%>%
  slice_max(n=15, order_by=sum)
# Phylum                      sum
# <chr>                     <dbl>
# 1 " p__Proteobacteria"    0.702  
# 2 " p__Bacteroidota"      0.0794 
# 3 " p__Actinobacteriota"  0.0425 
# 4 " p__Planctomycetota"   0.0409 
# 5 " p__WPS-2"             0.0280 
# 6 " p__Verrucomicrobiota" 0.0255 
# 7 " p__Nitrospirota"      0.0221 
# 8 " p__Acidobacteriota"   0.0148 
# 9 " p__Cyanobacteria"     0.0144 
# 10 " p__Chloroflexi"       0.0138 
# 11 " p__Gemmatimonadota"   0.00584
# 12 " p__Myxococcota"       0.00405
# 13 " p__Patescibacteria"   0.00372
# 14 " p__Armatimonadota"    0.00339

## Same for the family
taxglom_biogeo_family <- tax_glom(physeq=merged_NOMIS_core_asv, taxrank=rank_names(merged_NOMIS_core_asv)[5], NArm=F)
tax_table_core_family <- tax_table(taxglom_biogeo_family)
## Take the relative abundance and group by mountain ranges
transf_ab_fam = transform_sample_counts(taxglom_biogeo_family, function(x) x / sum(x))
sample_merge_fam = merge_samples(transf_ab_fam, "Site_c")
trans_ra_fam = transform_sample_counts(sample_merge_fam, function(x) x / sum(x))

asv_core_fam <- otu_table(trans_ra_fam, taxa_are_rows=T)
melt_core_fam <- psmelt(asv_core_fam)
merge_taxo_core_fam <- merge(melt_core_fam, tax_table_core_family, by.x="OTU", by.y="row.names")

sumtot_family <- merge_taxo_core_fam  %>% group_by(Family) %>% summarize(sum=sum(Abundance/10))%>%
  slice_max(n=15, order_by=sum)


## Same for the genera
taxglom_biogeo_genera <- tax_glom(physeq=merged_NOMIS_core_asv, taxrank=rank_names(merged_NOMIS_core_asv)[6], NArm=F)
tax_table_core_genera <- tax_table(taxglom_biogeo_genera)
## Take the relative abundance and group by mountain ranges
transf_ab_gen = transform_sample_counts(taxglom_biogeo_genera, function(x) x / sum(x))
sample_merge_gen = merge_samples(transf_ab_gen, "Site_c")
trans_ra_gen = transform_sample_counts(sample_merge_gen, function(x) x / sum(x))

asv_core_gen <- otu_table(trans_ra_gen, taxa_are_rows=T)
melt_core_gen <- psmelt(asv_core_gen)
merge_taxo_core <- merge(melt_core_gen, tax_table_core_genera, by.x="OTU", by.y="row.names")

##Now we would like to know what are the most abundant phyla across the dataset.
sumtot_genera <- merge_taxo_core %>% group_by(Genus) %>% summarize(sum=sum(Abundance/10))%>%
filter(!(Genus %in% c(""," g__uncultured"))) %>%
  slice_max(n=15, order_by=sum)
# Genus                          sum
# <chr>                        <dbl>
# 1 " g__Polaromonas"           0.113 
# 2 " g__Methylotenera"         0.105 
# 3 " g__Rhodoferax"            0.0671
# 4 " g__Sphingomonas"          0.0614
# 5 " g__Thiobacillus"          0.0438
# 6 " g__Ferruginibacter"       0.0365
# 7 " g__Rhizobacter"           0.0353
# 8 " g__Novosphingobium"       0.0334
# 9 " g__WPS-2"                 0.0280
# 10 " g__Ellin6067"             0.0257
# 11 " g__Parablastomonas"       0.0244
# 12 " g__Nitrospira"            0.0221
# 13 " g__Pirellula"             0.0168
# 14 " g__Arenimonas"            0.0154
# 15 " g__CL500-29_marine_group" 0.0153


 

