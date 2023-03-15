### Endemic, core and Character ASVs 

metadata_nomis <- as.data.frame(sample_data(prune_Uganda))
nomis_asv_count <- otu_table(prune_Uganda, taxa_are_rows=T)
nomis_asv_count_df <- nomis_asv_count [rowSums(nomis_asv_count [])>0,]

##Core taxa - count table
binary_tax_merge <- load("binary_tax_merge.RData")
merge_core_abondance <- merge(nomis_asv_count_df, binary_tax_merge, by.x="row.names", by.y="ASV" )
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

##Endemic taxa - count table
merge_endemic_phylo <- readRDS("merge_endemic_phylo.RDS")
asv_endemic <- otu_table(merge_endemic_phylo, taxa_are_rows=T)

##Character ASVs - count table 
character_asv <- read.csv("character_ASV_unique_20221412.csv")
merge_character_abondance <- merge(nomis_asv_count_df, character_asv, by.x="row.names", by.y="Q.ASVs" )
rownames(merge_character_abondance) <- merge_character_abondance$Row.names
merge_character_abondance$Row.names <- NULL
character_table <- otu_table(merge_character_abondance, taxa_are_rows=T)
merge_character_phylo <- merge_phyloseq(character_table, tax_core, metadata_glaciers)
saveRDS(merge_character_phylo, "merge_character_phylo_20221412.RDS")

##merging the 3 datasets// 
##removing what was redundant

unique_asv <- read.csv("unique_endemic_core_character_20221412.csv",sep=",",header=T)
nomis_asv_count_unique <- nomis_asv_count_df[rownames(nomis_asv_count_df) %in% unique_asv$ASV,]
asv_nomis_unique <- otu_table(nomis_asv_count_unique, taxa_are_rows=T)
nomis_asv_unique_phylo <- merge_phyloseq(asv_nomis_unique, tax_core, metadata_glaciers)

unique_Genus_taxglom <- tax_glom(nomis_asv_unique_phylo, taxrank=rank_names(nomis_asv_unique_phylo)[6], NArm=F)
transf_unique = transform_sample_counts(unique_Genus_taxglom, function(x) x / sum(x))

TopASV_f <- names(sort(taxa_sums(transf_unique), TRUE)[1:11])
top25_NOMIS_f <- prune_species(TopASV_f, transf_unique)
top25_NOMIS_f <- prune_taxa(taxa_sums(top25_NOMIS_f)>0, top25_NOMIS_f)
top_Genus<-as.data.frame(tax_table(top25_NOMIS_f))

merged_NOMIS_core_ab <- readRDS("merge_NOMIS_core_ab.RDS")
core_RA = transform_sample_counts(merged_NOMIS_core_ab, function(x) x / sum(x))
core_RA = merge_samples(core_RA, "Site_c")
core_RA = transform_sample_counts(core_RA, function(x) x / sum(x))

data_core <- psmelt(core_RA) # create dataframe from phyloseq object

data_core$Genus <-as.character(data_core$Genus) #convert to character

##We kept the 15 most abundant genera
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

##But lets find a nice color palette
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


## Same with endemic 
merge_endemic_phylo <- readRDS("merge_endemic_phylo.RDS")
endemic_RA = transform_sample_counts(merge_endemic_phylo, function(x) x / sum(x))
endemic_RA = merge_samples(endemic_RA, "Site_c")
endemic_RA = transform_sample_counts(endemic_RA, function(x) x / sum(x))

data_endemic <- psmelt(endemic_RA) # create dataframe from phyloseq object

data_endemic$Genus <-as.character(data_endemic$Genus) #convert to character

##We kept the 15 most abundant genera
sumtot_endemic <-
  data_endemic %>% group_by(Genus) %>% summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% top_Genus$Genus) %>% filter(!(Genus %in% c(""," g__uncultured")))
  
##Here we set to "other" 
data_endemic$Genus[!(data_endemic$Genus %in% sumtot_endemic$Genus)] <- "Other"

data_endemic$endemic <- "endemic"
data_endemic$totalAbundance <- sum(data_endemic$Abundance)

data_endemic_mod <- data_endemic%>%
  group_by(Genus, endemic)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_endemic <- ggplot(data=data_endemic_mod, aes(x=endemic, y=abundance, fill=Genus))
barplot_biogeo_endemic <- barplot_biogeo_endemic + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))

barplot_biogeo_endemic<- barplot_biogeo_endemic + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


## With Character ASVs
merge_character_phylo <- readRDS("merge_character_phylo_20221412.RDS")
character_RA = transform_sample_counts(merge_character_phylo, function(x) x / sum(x))
character_RA = merge_samples(character_RA, "Site_c")
character_RA = transform_sample_counts(character_RA, function(x) x / sum(x))

data_character <- psmelt(character_RA) # create dataframe from phyloseq object

data_character$Genus <-as.character(data_character$Genus) #convert to character

##We kept the 15 most abundant genera
sumtot_character <-
  data_character %>% group_by(Genus) %>% summarize(sum = sum(Abundance)) %>%
  filter(Genus %in% top_Genus$Genus) %>%filter(!(Genus %in% c(""," g__uncultured")))

##Here we set to "other" 
data_character$Genus[!(data_character$Genus %in% sumtot_character$Genus)] <- "Other"

data_character$character <- "character"
data_character$totalAbundance <- sum(data_character$Abundance)

data_character_mod <- data_character%>%
  group_by(Genus, character)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()


n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_character <- ggplot(data=data_character_mod, aes(x=character, y=abundance, fill=Genus))
barplot_biogeo_character <- barplot_biogeo_character + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
barplot_biogeo_character <- barplot_biogeo_character + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))

##Merge the 3 graphs together...with ggarange

ggarrange(barplot_biogeo_endemic,barplot_biogeo_core,barplot_biogeo_character, ncol = 3, nrow = 1, common.legend=T)

##ggrdiges - density plot

##D'abord on enlève les Abundance qui sont nulles pour tous les datasets
dataset_endemic_filter <- as.data.frame(data_endemic[data_endemic$Abundance >0,])
rename_endemic <- rename(dataset_endemic_filter, category = endemic)
dataset_core_filter <- as.data.frame(data_core[data_core$Abundance >0,])
dataset_core_filter <- as.data.frame(subset(dataset_core_filter, select=-c(X, X.1, X.2, X.3, X.4,X.5,X.6,X.7)))
rename_core <- rename(dataset_core_filter, category = core)
dataset_character_filter <- as.data.frame(data_character[data_character$Abundance >0,])
rename_character <- rename(dataset_character_filter, category = character)

binddataset <- rbind(rename_endemic, rename_core, rename_character)


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

##Here we are going to filter the taxa that are endemic and core from the character
##On merge les trois jeux de données puis on prendre l'abondance relative et voit quels sont
##les familles les plus abondantes.

##So we need to upload the core, character, and endemic list 
##Here are the metadata and the full ASV table
metadata_nomis <- as.data.frame(sample_data(prune_Uganda))
nomis_asv_count <- otu_table(prune_Uganda, taxa_are_rows=T)
nomis_asv_count_df <- nomis_asv_count [rowSums(nomis_asv_count [])>0,]

##Core taxa - count table
binary_tax_merge <- load("binary_tax_merge.RData")
merge_core_abondance <- merge(nomis_asv_count_df, binary_tax_merge, by.x="row.names", by.y="ASV" )
##We keep only the community matrix and remove everything else.
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

##Endemic taxa - count table
merge_endemic_phylo <- readRDS("merge_endemic_phylo.RDS")
asv_endemic <- otu_table(merge_endemic_phylo, taxa_are_rows=T)

##Character ASVs - count table 
character_asv <- read.csv("character_ASV_unique_20221412.csv")
merge_character_abondance <- merge(nomis_asv_count_df, character_asv, by.x="row.names", by.y="Q.ASVs" )
rownames(merge_character_abondance) <- merge_character_abondance$Row.names
merge_character_abondance$Row.names <- NULL
character_table <- otu_table(merge_character_abondance, taxa_are_rows=T)
merge_character_phylo <- merge_phyloseq(character_table, tax_core, metadata_glaciers)
saveRDS(merge_character_phylo, "merge_character_phylo_20221412.RDS")
##Mtn il faut merger les trois jeux de données // 
##On copie colle dans un fichier excel les ASVs -- total 34358 entries and we removed what was redundant
##Remaining 31857

unique_asv <- read.csv("unique_endemic_core_character_20221412.csv",sep=",",header=T)
nomis_asv_count_unique <- nomis_asv_count_df[rownames(nomis_asv_count_df) %in% unique_asv$ASV,]
asv_nomis_unique <- otu_table(nomis_asv_count_unique, taxa_are_rows=T)
nomis_asv_unique_phylo <- merge_phyloseq(asv_nomis_unique, tax_core, metadata_glaciers)

unique_Family_taxglom <- tax_glom(nomis_asv_unique_phylo, taxrank=rank_names(nomis_asv_unique_phylo)[5], NArm=F)
transf_unique = transform_sample_counts(unique_Family_taxglom, function(x) x / sum(x))

TopASV_f <- names(sort(taxa_sums(transf_unique), TRUE)[1:11])
top25_NOMIS_f <- prune_species(TopASV_f, transf_unique)
top25_NOMIS_f <- prune_taxa(taxa_sums(top25_NOMIS_f)>0, top25_NOMIS_f)
top_Family<-as.data.frame(tax_table(top25_NOMIS_f))

##Il nous faut la table du core sous forme d'abondance relative
merged_NOMIS_core_ab <- readRDS("merge_NOMIS_core_ab.RDS")
core_RA = transform_sample_counts(merged_NOMIS_core_ab, function(x) x / sum(x))
core_RA = merge_samples(core_RA, "Site_c")
core_RA = transform_sample_counts(core_RA, function(x) x / sum(x))

data_core <- psmelt(core_RA) # create dataframe from phyloseq object

data_core$Family <-as.character(data_core$Family) #convert to character

##We kept the 15 most abundant genera

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

##But lets find a nice color palette
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


## Same with endemic 
merge_endemic_phylo <- readRDS("merge_endemic_phylo.RDS")
endemic_RA = transform_sample_counts(merge_endemic_phylo, function(x) x / sum(x))
endemic_RA = merge_samples(endemic_RA, "Site_c")
endemic_RA = transform_sample_counts(endemic_RA, function(x) x / sum(x))

data_endemic <- psmelt(endemic_RA) # create dataframe from phyloseq object

data_endemic$Family <-as.character(data_endemic$Family) #convert to character

##We kept the 15 most abundant genera
sumtot_endemic <-
  data_endemic %>% group_by(Family) %>% summarize(sum = sum(Abundance)) %>%
  filter(Family %in% top_Family$Family) %>% filter(!(Family %in% c(""," g__uncultured")))

##Here we set to "other" 
data_endemic$Family[!(data_endemic$Family %in% sumtot_endemic$Family)] <- "Other"

data_endemic$endemic <- "endemic"
data_endemic$totalAbundance <- sum(data_endemic$Abundance)

data_endemic_mod <- data_endemic%>%
  group_by(Family, endemic)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_endemic <- ggplot(data=data_endemic_mod, aes(x=endemic, y=abundance, fill=Family))
barplot_biogeo_endemic <- barplot_biogeo_endemic + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))

barplot_biogeo_endemic<- barplot_biogeo_endemic + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                        panel.background = element_blank(), axis.line = element_line(colour = "black"))


## With Character ASVs
merge_character_phylo <- readRDS("merge_character_phylo_20221412.RDS")
character_RA = transform_sample_counts(merge_character_phylo, function(x) x / sum(x))
character_RA = merge_samples(character_RA, "Site_c")
character_RA = transform_sample_counts(character_RA, function(x) x / sum(x))

data_character <- psmelt(character_RA) # create dataframe from phyloseq object

data_character$Family <-as.character(data_character$Family) #convert to character

##We kept the 15 most abundant genera
sumtot_character <-
  data_character %>% group_by(Family) %>% summarize(sum = sum(Abundance)) %>%
  filter(Family %in% top_Family$Family) %>%filter(!(Family %in% c(""," g__uncultured")))

##Here we set to "other" 
data_character$Family[!(data_character$Family %in% sumtot_character$Family)] <- "Other"

data_character$character <- "character"
data_character$totalAbundance <- sum(data_character$Abundance)

data_character_mod <- data_character%>%
  group_by(Family, character)%>%
  summarise(abundance = sum(Abundance)/totalAbundance)%>%
  distinct()

n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo_character <- ggplot(data=data_character_mod, aes(x=character, y=abundance, fill=Family))
barplot_biogeo_character <- barplot_biogeo_character + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#4E79A7FF", "#A0CBE8FF", "#F28E2BFF", "#FFBE7DFF", "#59A14FFF", "#8CD17DFF", "#B6992DFF", 
                               "#F1CE63FF" ,"#499894FF", "#86BCB6FF", "#E15759FF", "#FF9D9AFF", "#79706EFF", "#BAB0ACFF","#B07AA1FF", "#D4A6C8FF", "#9D7660FF", "#D7B5A6FF" )) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
barplot_biogeo_character <- barplot_biogeo_character + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))

##Merge the 3 graphs together...with ggarange

ggarrange(barplot_biogeo_endemic,barplot_biogeo_core,barplot_biogeo_character, ncol = 3, nrow = 1, common.legend=T)

##ggrdiges - density plot
dataset_endemic_filter <- as.data.frame(data_endemic[data_endemic$Abundance >0,])
rename_endemic <- rename(dataset_endemic_filter, category = endemic)
dataset_core_filter <- as.data.frame(data_core[data_core$Abundance >0,])
dataset_core_filter <- as.data.frame(subset(dataset_core_filter, select=-c(X, X.1, X.2, X.3, X.4,X.5,X.6,X.7)))
rename_core <- rename(dataset_core_filter, category = core)
dataset_character_filter <- as.data.frame(data_character[data_character$Abundance >0,])
rename_character <- rename(dataset_character_filter, category = character)

binddataset <- rbind(rename_endemic, rename_core, rename_character)


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

