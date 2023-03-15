### Loading NOMIS data 
library(phyloseq)
library(vegan)
library(data.table)
library(speedyseq)
library(phyloseqCompanion)
library(dplyr)
asv<-fread(file="1022_NOMIS_16S_merged_table.csv",header=TRUE, sep=",")
newnames <- colnames(asv)
#\\d -> nombre et $:la fin - ce que l'on cherche ce sont toutes les valeurs qui finissent par un nombre. les parenthèses c'est pour dire que ce qui est entre parenthèses
#est gardé -- le 1 -> c'est pour le premier set de parenthèse
newnames <- gsub("(.*_\\d)$","\\1_16S_sed",newnames)
asv_df<-as.data.frame(asv)
colnames(asv_df) <- newnames
rownames(asv_df) <- asv_df$Feature_ID
asv_df$Feature_ID <- NULL

tax<-read.csv(file="1022_NOMIS_16S_merged_taxonomy.csv",sep=",",header=TRUE,row.names=1)

## Import QIIME pipeline map-file---ALWAYS ALWAYS put int in txt. Uganda samples are now included
metadata_NOMIS="metadata_NOMIS_1122_FULL.txt"
metadata_NOMIS<-import_qiime_sample_data(metadata_NOMIS)

OTU_NOMIS <- otu_table(asv_df,taxa_are_rows=TRUE)
tax_NOMIS <- tax_table(as.matrix(tax))

## Merge everything (mapping file, tree and OTUs)
merged_NOMIS <- merge_phyloseq(OTU_NOMIS, tax_NOMIS, metadata_NOMIS)

nomis_meta <- sample.data.frame(merged_NOMIS)

## Subsetting taxa 
merged_NOMIS_sub <- subset_taxa(merged_NOMIS, (Kingdom!="d__Eukaryota") | is.na(Kingdom)) 
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_sub, (Kingdom!="d__Archaea") | is.na(Kingdom))
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_sub, (Order!=" o__Chloroplast") )
merged_NOMIS_sub <- subset_taxa(merged_NOMIS_sub, (Family!=" f__Mitochondria"))

## prune doubletons
pruned_doubletons <- prune_taxa(taxa_sums(merged_NOMIS_sub) > 2, merged_NOMIS_sub) 

## Select Substrate and Position
site = c("sed")
sediment <- subset_samples(pruned_doubletons,substrate %in% site)

posi = c("U")
sediment_up <- subset_samples(sediment, Position %in% posi)

## Rarefaction curves
set.seed(42)
calculate_rarefaction_curves <- function(sediment_up, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(sediment_up, measures, depth) {
    if(max(sample_sums(sediment_up)) < depth) return()
    sediment_up <- prune_samples(sample_sums(sediment_up) >= depth, sediment_up)
    
    rarified_sediment_up <- rarefy_even_depth(sediment_up, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_sediment_up, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, sediment_up = sediment_up, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(sediment_up, c('Observed', 'Shannon'), rep(c(1, 10, 100, 1000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))

rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(sediment_up)), by.x = 'Sample', by.y = 'row.names')

rarefaction_curve_observed <- subset(rarefaction_curve_data_summary_verbose, Measure == "Observed")
rarefaction_curve_shannon <- subset(rarefaction_curve_data_summary_verbose, Measure == "shannon")


library('ggplot2')
raref_plot <- ggplot(
  data = rarefaction_curve_observed,
  mapping = aes(
    x = Depth,
    y = Alpha_diversity_mean,
    ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd,
    colour =Site_c,
    group = Sample
  )
) + geom_line(size=0.2
) + facet_wrap(
  facets = ~ Site_c
  )

raref_plot + scale_color_manual(values= c("Alaska"="#2E2A2BFF","Alps"="#CF4E9CFF","Caucasus"="#8C57A2FF",
                                              "Chile"="#3EBCB6","Ecuador"="#82581FFF","Greenland"="#2F509EFF",
                                              "kyrgyzstan"="#E5614CFF","Nepal"="#97A1A7FF","New_Zealand"="#bee183","Norway"="#DC9445FF","Uganda"="#EDD03E"))+ theme_bw() +  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

+ geom_pointrange(size=0.5
## Filtering taxa for prevalence 
## Taxa should be present in at least 2 out of 3 biological replicates or 1/2. Creating a table of ASVs
set.seed(132)
filtering_br <- function(phyloseq_obj){
  meta_sediment <- as.data.frame(sample_data(phyloseq_obj))
  asv_sedimentup <- as.data.frame(otu_table(phyloseq_obj, taxa_are_rows = T))
  
  new_asv <- c()
  for (Sample in meta_sediment$code_gl) {
    sub_asvsed <- asv_sedimentup[,startsWith(colnames(asv_sedimentup), paste0(Sample,'_',collapse = ''))]
    if (is.null(sub_asvsed)){}
    else{
      if (!(is.vector(sub_asvsed))){
        sub_asv_preval <- rowMeans(sub_asvsed > 0)
        new_asv = c(new_asv, rownames(sub_asvsed)[sub_asv_preval >= 0.5])}
      else {
        new_asv = c(new_asv, names(sub_asvsed)[sub_asvsed > 0])}
      }}
  new_asv = unique(new_asv)
 return(new_asv)
}
##List of ASVs
asvsubset<-filtering_br(sediment_up)
saveRDS(asvsubset, "asvsubset.RData")

##Now we subset the ASVs filtered from the ASV table
asvfull <- as.data.frame(otu_table(sediment_up, taxa_are_rows=T))
subasv <- asvfull[rownames(asvfull) %in% asvsubset,]
subasv_table <- otu_table(subasv, taxa_are_rows=T)
saveRDS(subasv_table, "subasv_table.RData")
## Melt ASV_table
melt_asv <- melt(subasv_table)

code_glaciers <- sample.data.frame(metadata_NOMIS[,c("Sample","code_gl")])
rownames(code_glaciers)<-NULL

## Merge ASV_table and code_gl
melt_merge_as <- merge(melt_asv, code_glaciers, by.x="Var2", by.y="Sample")
saveRDS(melt_merge_as,"20221026melt_merge.RData")

melt_merge_asv <- melt_merge_as %>%
  group_by(code_gl, Var1) %>% 
  summarise(average=mean(value)) %>% mutate(across(where(is.numeric), ~ as.integer(round(.x))))

dcast_asv <- dcast(melt_merge_asv, formula= Var1 ~ code_gl)
rownames(dcast_asv) <- dcast_asv$Var1
dcast_asv$Var1 <- NULL
asv_table_dcast <- as.matrix(otu_table(dcast_asv, taxa_are_rows=T))

metadata_glaciers="metadata_glaciers_nomis_2022612.tsv"
metadata_glaciers<-import_qiime_sample_data(metadata_glaciers)

tax_sed <- as.matrix(tax_table(sediment_up))

merge_glaciers_full <- merge_phyloseq(asv_table_dcast, tax_sed, metadata_glaciers)
unrarefied_asv <- otu_table(merge_glaciers_full, taxa_are_rows=T)
unrarefied_taxtable <- tax_table(merge_glaciers_full)
write.csv(unrarefied_asv, "2022114_unrarefied_asv.csv")
write.csv(unrarefied_taxtable, "2022114_unrarefied_taxtable.csv")

##Rarefaction
NOMIS_FR <- rarefy_even_depth(merge_glaciers_full, sample.size=min(sample_sums(merge_glaciers_full)), rngseed=678, replace=F, trimOTUs=TRUE, verbose=TRUE)
#NOMIS_FR <- saveRDS(NOMIS_FR,"NOMIS_FR_20230901.RData")
NOMIS_FR <- readRDS("NOMIS_FR_20230901.RData")

NOMIS_full_asv_sedimentup <- otu_table(NOMIS_FR, taxa_are_rows=T)
NOMIS_full_taxonomy_sedimentup <- tax_table(NOMIS_FR)
#saveRDS(NOMIS_FR,"NOMIS_FR.RData")

tax_for_tree<-as.data.frame(tax_table(NOMIS_FR))
first_row <- rownames(tax_for_tree)
write.csv(first_row,"20221026_tree_ASV_NOMIS.csv")
write.csv(NOMIS_full_asv_sedimentup,"20221026NOMIS_full_asv_sedimentup_filtered_final.csv")
write.csv(NOMIS_full_taxonomy_sedimentup,"20221026NOMIS_full_tax_sedimentup_filtered_final.csv")







