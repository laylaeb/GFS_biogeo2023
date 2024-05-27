library(speedyseq)
library(phyloseq)
library(phyloseqCompanion)

# Specific, core and indicator taxa - relative abundance
# Need to upload the 3 lists of taxa
# Some metrics regarding the Number of ASVs found across the dataset//Removing the Uganda sample!

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

## remove Uganda because only 1 GFS
uganda=c("Uganda")
prune_Uganda <- subset_samples(merged_NOMIS_DEBLUR , !site_c %in% uganda)
nomis_asv_count_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))
nomis_asv_count_df <- nomis_asv_count_df[rowSums(nomis_asv_count_df)> 0, ] ### 54019 ASVs and 151 GFSs

## Core taxa - count table
core_list <-cores <- read_csv("20240301_core.csv")
merge_core_abondance <- nomis_asv_count_df[rownames(nomis_asv_count_df) %in% core_list$ASV,]

mca_table <- otu_table(merge_core_abondance, taxa_are_rows=T)

## Merge everything -- CORE ASV
merged_NOMIS_core_ab<- merge_phyloseq(mca_table, tax_NOMIS, metadata_NOMIS)
#saveRDS(merged_NOMIS_core_ab, "2024_merge_NOMIS_core_ab.RDS")

# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 165 taxa and 151 samples ]:
# sample_data() Sample Data:        [ 151 samples by 5 sample variables ]:
# tax_table()   Taxonomy Table:     [ 165 taxa by 7 taxonomic ranks ]:
# taxa are rows

## Specific taxa - count table
merge_Specific <- read_csv("202403_Specific_list.csv")

merge_Specific_abundance <- nomis_asv_count_df[rownames(nomis_asv_count_df) %in% merge_Specific$ASV,]

asv_Specific_table <- otu_table(merge_Specific_abundance, taxa_are_rows=T)
merge_NOMIS_Specific_ab <- merge_phyloseq(asv_Specific_table, tax_NOMIS, metadata_NOMIS)
#saveRDS(merge_NOMIS_Specific_ab, "2024_merge_NOMIS_Specific_ab.RDS")

## Indicator ASVs
indicator <- read_csv("2024_0307_indicator_taxa.csv")
indicator$ASV <- indicator$...1
indicator$...1 <- NULL

merge_indicator_abundance <- nomis_asv_count_df[rownames(nomis_asv_count_df) %in% indicator$ASV,]

indicator_table <- otu_table(merge_indicator_abundance, taxa_are_rows=T)

## We need to filter ASVs that are not present in core and Specific taxa 
indicator_table_filtered <- subset(indicator_table, !(row.names(indicator_table) %in% c(row.names(asv_Specific_table), row.names(mca_table))))
merge_indicator_phylo <- merge_phyloseq(indicator_table_filtered, tax_NOMIS, metadata_NOMIS)

#saveRDS(merge_indicator_phylo, "2024_merge_indicator_phylo.RDS")
#indicator_table_filtered <- readRDS("2024_merge_indicator_phylo.RDS")

indicator_ASV <- row.names(indicator_table_filtered)
Specific_ASV <- row.names(asv_Specific_table)
core_ASV <- row.names(mca_table)

df_core_Specific_indicator <- c(indicator_ASV,Specific_ASV,core_ASV)
##40418 elements merged!
df_core_Specific_indicator_unique <- unique(df_core_Specific_indicator)
##40418 elements

## Now we do have 40418 ASVs we will subset them to the full ASV table and plot the taxonomy.
nomis_asv_count_unique <- nomis_asv_count_df[rownames(nomis_asv_count_df) %in% df_core_Specific_indicator_unique,]
asv_nomis_unique <- otu_table(nomis_asv_count_unique, taxa_are_rows=T)
nomis_asv_unique_phylo <- merge_phyloseq(asv_nomis_unique, tax_NOMIS, metadata_NOMIS)

unique_Genus_taxglom <- tax_glom(nomis_asv_unique_phylo, taxrank=rank_names(nomis_asv_unique_phylo)[6], NArm=F)
transf_unique = transform_sample_counts(unique_Genus_taxglom, function(x) x / sum(x))

unique_phylum_taxglom <- tax_glom(nomis_asv_unique_phylo, taxrank=rank_names(nomis_asv_unique_phylo)[2], NArm=F)
transf_unique_phylum = transform_sample_counts(unique_phylum_taxglom, function(x) x / sum(x))

unique_family_taxglom <- tax_glom(nomis_asv_unique_phylo, taxrank=rank_names(nomis_asv_unique_phylo)[5], NArm=F)
transf_family = transform_sample_counts(unique_family_taxglom, function(x) x / sum(x))

TopASV_f <- names(sort(taxa_sums(transf_family), TRUE)[1:11])
top25_NOMIS_f <- prune_species(TopASV_f, transf_unique)
top25_NOMIS_f <- prune_taxa(taxa_sums(top25_NOMIS_f)>0, top25_NOMIS_f)
top_Genus<-as.data.frame(tax_table(top25_NOMIS_f))

asv_table_all_genus <- otu_table(transf_unique, taxa_are_rows=T)
tax_table_all_genus <- tax_table(transf_unique)

asv_table_all_family <- otu_table(transf_family, taxa_are_rows=T)
tax_table_all_family <- tax_table(transf_family)

asv_table_all_phylum <- otu_table(transf_unique_phylum, taxa_are_rows=T)
tax_table_all_phylum <- tax_table(transf_unique_phylum)

write.csv(asv_table_all_genus, "asv_table_all_genus.csv")
write.csv(tax_table_all_genus, "tax_table_all_genus.csv")
# 
write.csv(asv_table_all_family, "asv_table_all_family.csv")
write.csv(tax_table_all_family, "tax_table_all_family.csv")
# 
write.csv(asv_table_all_phylum, "asv_table_all_phylum.csv")
write.csv(tax_table_all_phylum, "tax_table_all_phylum.csv")

## We do have a list of the genera to plot
## Lets start with the core microbiome

## computing core - relative abundance
merged_NOMIS_core_ab <- readRDS("merge_NOMIS_core_ab.RDS")
core_RA = transform_sample_counts(merged_NOMIS_core_ab, function(x) x / sum(x))
core_RA = merge_samples(core_RA, "site_c")
core_RA = transform_sample_counts(core_RA, function(x) x / sum(x))

data_core <- psmelt(core_RA) # create dataframe from phyloseq object
data_core$Genus <-as.character(data_core$Genus) # convert to character

## We kept the 15 most abundant genera
sumtot_core <-
  data_core %>% group_by(Genus) %>% summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% top_Genus$Genus) %>%filter(!(Genus %in% c(""," g__uncultured")))

data_core$Genus[!(data_core$Genus %in% sumtot_core$Genus)] <- "Other"

data_core$core <- "core"
data_core$totalAbundance <- sum(data_core$Abundance)

data_core_mod <- data_core%>%
  group_by(Genus, core)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_core <- ggplot(data=data_core_mod, aes(x=core, y=abundance, fill=Genus))
barplot_biogeo_core <- barplot_biogeo_core + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
barplot_biogeo_core<- barplot_biogeo_core+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))

## Let's do the same with the Specific ASVs
merge_Specific_phylo <- readRDS("merge_Specific_phylo.RDS")
Specific_RA = transform_sample_counts(merge_Specific_phylo, function(x) x / sum(x))
Specific_RA = merge_samples(Specific_RA, "site_c")
Specific_RA = transform_sample_counts(Specific_RA, function(x) x / sum(x))

data_Specific <- psmelt(Specific_RA) # create dataframe from phyloseq object

data_Specific$Genus <-as.character(data_Specific$Genus) #convert to character

## We kept the 15 most abundant genera
sumtot_Specific <-
  data_Specific %>% group_by(Genus) %>% summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% top_Genus$Genus) %>% filter(!(Genus %in% c(""," g__uncultured")))
  
## Here we set to "other" 
data_Specific$Genus[!(data_Specific$Genus %in% sumtot_Specific$Genus)] <- "Other"

data_Specific$Specific <- "Specific"
data_Specific$totalAbundance <- sum(data_Specific$Abundance)

data_Specific_mod <- data_Specific%>%
  group_by(Genus, Specific)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_Specific <- ggplot(data=data_Specific_mod, aes(x=Specific, y=abundance, fill=Genus))
barplot_biogeo_Specific <- barplot_biogeo_Specific + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))

barplot_biogeo_Specific<- barplot_biogeo_Specific + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


## With indicator ASVs
merge_indicator_phylo <- readRDS("merge_indicator_phylo_2023.RDS")
indicator_RA = transform_sample_counts(merge_indicator_phylo, function(x) x / sum(x))
indicator_RA = merge_samples(indicator_RA, "site_c")
indicator_RA = transform_sample_counts(indicator_RA, function(x) x / sum(x))

data_indicator <- psmelt(indicator_RA) # create dataframe from phyloseq object

data_indicator$Genus <-as.character(data_indicator$Genus) #convert to indicator

## We kept the 15 most abundant genera
sumtot_indicator <-
  data_indicator %>% group_by(Genus) %>% summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% top_Genus$Genus) %>%filter(!(Genus %in% c(""," g__uncultured")))

## Here we set to "other" 
data_indicator$Genus[!(data_indicator$Genus %in% sumtot_indicator$Genus)] <- "Other"

data_indicator$indicator <- "indicator"
data_indicator$totalAbundance <- sum(data_indicator$Abundance)

data_indicator_mod <- data_indicator%>%
  group_by(Genus, indicator)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()


n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_indicator <- ggplot(data=data_indicator_mod, aes(x=indicator, y=abundance, fill=Genus))
barplot_biogeo_indicator <- barplot_biogeo_indicator + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
barplot_biogeo_indicator <- barplot_biogeo_indicator + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))

## Merge the 3 graphs 
ggarrange(barplot_biogeo_Specific,barplot_biogeo_core,barplot_biogeo_indicator, ncol = 3, nrow = 1, common.legend=T)

## ggrdiges - density plot
dataset_Specific_filter <- as.data.frame(data_Specific[data_Specific$Abundance >0,])
rename_Specific <- rename(dataset_Specific_filter, category = Specific)

dataset_core_filter <- as.data.frame(data_core[data_core$Abundance >0,])
rename_core <- rename(dataset_core_filter, category = core)

dataset_indicator_filter <- as.data.frame(data_indicator[data_indicator$Abundance >0,])
rename_indicator <- rename(dataset_indicator_filter, category = indicator)

binddataset <- rbind(rename_Specific, rename_core, rename_indicator)

##3 layers
library(tidyverse)
library(ggridges)
set.seed(3467)
# 
# binddataset_v2 <- binddataset %>%
#   mutate(Genus = fct_relevel(Genus, levels = " g__Ellin6067", " g__Ferruginibacter", " g__Flavobacterium"," g__Methylotenera"," g__Novosphingobium"," g__Polaromonas"," g__Rhizobacter"," g__Rhodoferax"," g__Sphingomonas"," g__Thiobacillus","Other"))
# 
# genus <- c(" g__Ellin6067", " g__Ferruginibacter", " g__Flavobacterium"," g__Methylotenera"," g__Novosphingobium"," g__Polaromonas"," g__Rhizobacter"," g__Rhodoferax"," g__Sphingomonas"," g__Thiobacillus","Other")
# 
# binddataset$Genus <- factor(binddataset$Genus)
# fct_reorder(binddataset$Genus, binddataset$Abundance, min)


ggplot(binddataset, aes(x = Abundance, y = fct_reorder(Genus, Abundance, .desc = F), fill=category)) + 
  geom_density_ridges(scale = 1, alpha=0.8) + 
  #scale_fill_cyclical(values = c("blue", "green","red"))+
  theme_bw() + scale_x_continuous(trans="log10") + scale_fill_brewer(palette = "Dark2")

############### Family

## Here we are going to filter the taxa that are Specific and core from the indicator
## We merge the dataset together and take the relative abundance
## Most abundant families 

## So we need to upload the core, indicator,and specific lists 
## Here are the metadata and the full ASV table
metadata_nomis <- as.data.frame(sample_data(prune_Uganda))
nomis_asv_count <- otu_table(prune_Uganda, taxa_are_rows=T)
nomis_asv_count_df <- nomis_asv_count [rowSums(nomis_asv_count [])>0,]

## Core taxa - count table
binary_tax_merge <- load("/Users/ezzat/Google Drive/Projet SBER/NOMIS/NOMIS_16S/NOMIS_October2022/binary_tax_merge.RData")
merge_core_abondance <- merge(nomis_asv_count_df, binary_tax_merge, by.x="row.names", by.y="ASV" )
## We keep only the community matrix and remove everything else.
write.csv(merge_core_abondance, "merge_core_abundance.csv")
merge_core_abondance <- read.csv("merge_core_abundance.csv",sep=",",header=T)
rownames(merge_core_abondance) <- merge_core_abondance$X
merge_core_abondance$X <- NULL
merge_core_ab <- otu_table(merge_core_abondance, taxa_are_rows=T)
tax_NOMIS <- tax_table(tax_core)
## Merge everything -- CORE ASV
merged_NOMIS_core_ab<- merge_phyloseq(merge_core_ab, tax_NOMIS, metadata_nomis)
saveRDS(merged_NOMIS_core_ab, "merge_NOMIS_core_ab.RDS")
# phyloseq-class experiment-level object
# otu_table()   OTU Table:          [ 200 taxa and 147 samples ]:
#   sample_data() Sample Data:        [ 147 samples by 48 sample variables ]:
#   tax_table()   Taxonomy Table:     [ 200 taxa by 7 taxonomic ranks ]:
#   taxa are rows

## Specific ASVs - count table
merge_Specific_phylo <- readRDS("merge_Specific_phylo.RDS")
asv_Specific <- otu_table(merge_Specific_phylo, taxa_are_rows=T)

## indicator ASVs - count table 
indicator_asv <- read.csv("indicator_ASV_unique_20221412.csv")
merge_indicator_abondance <- merge(nomis_asv_count_df, indicator_asv, by.x="row.names", by.y="Q.ASVs" )
rownames(merge_indicator_abondance) <- merge_indicator_abondance$Row.names
merge_indicator_abondance$Row.names <- NULL
indicator_table <- otu_table(merge_indicator_abondance, taxa_are_rows=T)
merge_indicator_phylo <- merge_phyloseq(indicator_table, tax_core, metadata_glaciers)
saveRDS(merge_indicator_phylo, "merge_indicator_phylo_20221412.RDS")

## to merge the 3 datasets

unique_asv <- read.csv("unique_Specific_core_indicator_20221412.csv",sep=",",header=T)
## On subset ces taxa
nomis_asv_count_unique <- nomis_asv_count_df[rownames(nomis_asv_count_df) %in% unique_asv$ASV,]
asv_nomis_unique <- otu_table(nomis_asv_count_unique, taxa_are_rows=T)
nomis_asv_unique_phylo <- merge_phyloseq(asv_nomis_unique, tax_core, metadata_glaciers)

unique_Family_taxglom <- tax_glom(nomis_asv_unique_phylo, taxrank=rank_names(nomis_asv_unique_phylo)[5], NArm=F)
transf_unique = transform_sample_counts(unique_Family_taxglom, function(x) x / sum(x))

TopASV_f <- names(sort(taxa_sums(transf_unique), TRUE)[1:11])
top25_NOMIS_f <- prune_species(TopASV_f, transf_unique)
top25_NOMIS_f <- prune_taxa(taxa_sums(top25_NOMIS_f)>0, top25_NOMIS_f)
top_Family<-as.data.frame(tax_table(top25_NOMIS_f))

# Core microbiome
# Need the table of the core with the relative abundance
merged_NOMIS_core_ab <- readRDS("merge_NOMIS_core_ab.RDS")
core_RA = transform_sample_counts(merged_NOMIS_core_ab, function(x) x / sum(x))
core_RA = merge_samples(core_RA, "Site_c")
core_RA = transform_sample_counts(core_RA, function(x) x / sum(x))

data_core <- psmelt(core_RA) # create dataframe from phyloseq object

data_core$Family <-as.indicator(data_core$Family) #convert to indicator

## We kept the 15 most abundant genera

sumtot_core <-
  data_core %>% group_by(Family) %>% summarize(sum = sum(Abundance)) %>%
  filter(Family %in% top_Family$Family) %>%filter(!(Family %in% c(""," g__uncultured")))

data_core$Family[!(data_core$Family %in% sumtot_core$Family)] <- "Other"

data_core$core <- "core"
data_core$totalAbundance <- sum(data_core$Abundance)

data_core_mod <- data_core%>%
  group_by(Family, core)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_core <- ggplot(data=data_core_mod, aes(x=core, y=abundance, fill=Family))
barplot_biogeo_core <- barplot_biogeo_core + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
barplot_biogeo_core<- barplot_biogeo_core+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))

## Same with Specific 
merge_Specific_phylo <- readRDS("merge_Specific_phylo.RDS")
Specific_RA = transform_sample_counts(merge_Specific_phylo, function(x) x / sum(x))
Specific_RA = merge_samples(Specific_RA, "Site_c")
Specific_RA = transform_sample_counts(Specific_RA, function(x) x / sum(x))

data_Specific <- psmelt(Specific_RA) # create dataframe from phyloseq object

data_Specific$Family <-as.indicator(data_Specific$Family) #convert to indicator

## We kept the 15 most abundant genera
sumtot_Specific <-
  data_Specific %>% group_by(Family) %>% summarize(sum = sum(Abundance)) %>%
  filter(Family %in% top_Family$Family) %>% filter(!(Family %in% c(""," g__uncultured")))

## Here we set to "other" 
data_Specific$Family[!(data_Specific$Family %in% sumtot_Specific$Family)] <- "Other"

data_Specific$Specific <- "Specific"
data_Specific$totalAbundance <- sum(data_Specific$Abundance)

data_Specific_mod <- data_Specific%>%
  group_by(Family, Specific)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_Specific <- ggplot(data=data_Specific_mod, aes(x=Specific, y=abundance, fill=Family))
barplot_biogeo_Specific <- barplot_biogeo_Specific + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))

barplot_biogeo_Specific<- barplot_biogeo_Specific + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


## With indicator ASVs
merge_indicator_phylo <- readRDS("merge_indicator_phylo_20221412.RDS")
indicator_RA = transform_sample_counts(merge_indicator_phylo, function(x) x / sum(x))
indicator_RA = merge_samples(indicator_RA, "Site_c")
indicator_RA = transform_sample_counts(indicator_RA, function(x) x / sum(x))

data_indicator <- psmelt(indicator_RA) # create dataframe from phyloseq object

data_indicator$Family <-as.indicator(data_indicator$Family) #convert to indicator

## We kept the 15 most abundant genera
sumtot_indicator <-
  data_indicator %>% group_by(Family) %>% summarize(sum = sum(Abundance)) %>%
  filter(Family %in% top_Family$Family) %>%filter(!(Family %in% c(""," g__uncultured")))

## Here we set to "other" 
data_indicator$Family[!(data_indicator$Family %in% sumtot_indicator$Family)] <- "Other"

data_indicator$indicator <- "indicator"
data_indicator$totalAbundance <- sum(data_indicator$Abundance)

data_indicator_mod <- data_indicator%>%
  group_by(Family, indicator)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_indicator <- ggplot(data=data_indicator_mod, aes(x=indicator, y=abundance, fill=Family))
barplot_biogeo_indicator <- barplot_biogeo_indicator + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
barplot_biogeo_indicator <- barplot_biogeo_indicator + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))

## Merge the 3 graphs together...with ggarange
ggarrange(barplot_biogeo_Specific,barplot_biogeo_core,barplot_biogeo_indicator, ncol = 3, nrow = 1, common.legend=T)

## ggrdiges - density plot
dataset_Specific_filter <- as.data.frame(data_Specific[data_Specific$Abundance >0,])
rename_Specific <- rename(dataset_Specific_filter, category = Specific)
dataset_core_filter <- as.data.frame(data_core[data_core$Abundance >0,])
dataset_core_filter <- as.data.frame(subset(dataset_core_filter, select=-c(X, X.1, X.2, X.3, X.4,X.5,X.6,X.7)))
rename_core <- rename(dataset_core_filter, category = core)
dataset_indicator_filter <- as.data.frame(data_indicator[data_indicator$Abundance >0,])
rename_indicator <- rename(dataset_indicator_filter, category = indicator)

binddataset <- rbind(rename_Specific, rename_core, rename_indicator)


##3 layers
library(tidyverse)
library(ggridges)
set.seed(3467)
# 
# binddataset_v2 <- binddataset %>%
#   mutate(Family = fct_relevel(Family, levels = " g__Ellin6067", " g__Ferruginibacter", " g__Flavobacterium"," g__Methylotenera"," g__Novosphingobium"," g__Polaromonas"," g__Rhizobacter"," g__Rhodoferax"," g__Sphingomonas"," g__Thiobacillus","Other"))
# 
# Family <- c(" g__Ellin6067", " g__Ferruginibacter", " g__Flavobacterium"," g__Methylotenera"," g__Novosphingobium"," g__Polaromonas"," g__Rhizobacter"," g__Rhodoferax"," g__Sphingomonas"," g__Thiobacillus","Other")
# 
# binddataset$Family <- factor(binddataset$Family)
# fct_reorder(binddataset$Family, binddataset$Abundance, min)


ggplot(binddataset, aes(x = Abundance, y = fct_reorder(Family, Abundance, .desc = F), fill =category)) + 
  geom_density_ridges(scale = 1, alpha=0.5) + 
  #scale_fill_cyclical(values = c("blue", "green","red"))+
  theme_bw() + scale_x_continuous(trans="log10") + scale_fill_brewer(palette = "Dark2")

