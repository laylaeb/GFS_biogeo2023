library(tidyverse)
library(purrr)
library(ape)
library(vegan)
library(ggpubr)
library(foreach)
library(doMC)
library(reshape2)

###############################################################
# Functions

filtering_br <- function(asv_tab, meta_tab){
  new_asv <- c()
  for (sample in meta_tab$sample) {
    sub_asvsed <- asv_tab[,colnames(asv_tab) == sample]
    if (!(is.vector(sub_asvsed))){
      sub_asv_preval <- rowMeans(sub_asvsed > 0)
      new_asv = c(new_asv, rownames(sub_asvsed)[sub_asv_preval >= 0.5])}
    else {
      new_asv = c(new_asv, names(sub_asvsed)[sub_asvsed > 0])}}
  new_asv = unique(new_asv)
  return(new_asv)}

rarefy_tab <- function(x, n){
  # Wrapper around the vegan rrarefy function (removes ASVs with no observations after rarefaction)
  rarefied_tab <- as.data.frame(as.data.frame(t(rrarefy(t(x), n))))
  rarefied_tab = rarefied_tab[rowSums(rarefied_tab) > 0,]
  return(rarefied_tab)}

take_avg <- function(x){
  # Function that creates an ASV table with summed replicates for each sample
  new_x = data.frame(row.names = rownames(x))
  # main for loop to create new data
  for (sample in unique(colnames(x))){
    sample_data = as.data.frame(x[,colnames(x) == sample])
    if (ncol(sample_data) > 1){new_x = cbind(new_x, rowMeans(sample_data))}
    else {new_x = cbind(new_x, sample_data)}}
  colnames(new_x) = unique(colnames(x))
  # remove ASVs with no observations after rarefaction
  new_x = new_x[rowSums(new_x) > 0,]
  # Set lowest count to 1
  #min_non_zero = min(new_x[new_x > 0])
  #new_x = min_non_zero * as.data.frame(new_x)
  #new_x = as.data.frame(sapply(new_x, as.integer))
  #melt_x = melt(as.matrix(x))
  #avg_x <- melt_x %>%
  #  group_by(Var1, Var2) %>% 
  #  summarise(average=mean(value)) %>% filter(average > 0)
  #dcast_x <- as.data.frame(dcast(avg_x, formula= Var1 ~ Var2, value.var = 'average', fill = 0))
  #
  #rownames(dcast_x) = dcast_x$Var1
  #dcast_x$Var1 = NULL
  
  min_non_zero = min(new_x[new_x > 0])
  new_x = as.data.frame(new_x) / min_non_zero
  new_x = as.data.frame(sapply(new_x, round))
  new_x = as.data.frame(sapply(new_x, as.integer))
  return(new_x)}

filter_taxa_samples <- function(tab, tax){
  # filter taxa
  if ('Taxon' %in% colnames(tax)){
    ASVs_to_keep = tax %>% filter(!grepl("d__Eukaryota", Taxon) , !grepl("d__;", Taxon), !grepl("d__Archaea", Taxon),
                                  !grepl(" o__Chloroplast", Taxon) , !grepl(" f__Mitochondria", Taxon)) %>% select(Feature_ID)}
  else{ASVs_to_keep = tax %>% filter(Kingdom != "d__Eukaryota", Kingdom != '', Kingdom != "d__Archaea", 
                                     Order != " o__Chloroplast", Family != " f__Mitochondria") %>% select(Feature_ID)}
  ASVs_to_remove = rownames(tab)[rowSums(tab[,colnames(tab) %in% c("NC_16S_KAUST2_rock","BLEX1","BLEX2","oldBL1","oldBL2")] > 0)]
  new_tab = tab[!(rownames(tab) %in% ASVs_to_remove),]
  new_tab = tab[rownames(tab) %in% ASVs_to_keep$Feature_ID,]
  
  # filter samples
  samples_to_keep = colnames(new_tab)[(grepl('_UP_', colnames(new_tab)) & grepl('16S_sed', colnames(new_tab)))]
  new_tab = new_tab[,colnames(new_tab) %in% samples_to_keep]
  
  # remove empty ASVs
  new_tab = new_tab[rowSums(new_tab) > 0,]
  return(new_tab)}

load_asv_data <- function(min_reads_cutoff){
  setwd('~/Desktop/leila_data')
  asv_tab_deblur = read.csv('data/081623_NOMIS_16S_merged_deblur_table.csv', sep=',')
  asv_tax_deblur = read.csv('data/081623_NOMIS_16S_merged_deblur_taxonomy.csv', sep=',')
  
  rownames(asv_tab_deblur) = asv_tab_deblur$Feature_ID
  
  asv_tab_deblur = filter_taxa_samples(asv_tab_deblur, asv_tax_deblur)
  asv_tab_deblur = asv_tab_deblur[,endsWith(colnames(asv_tab_deblur), 'sed') & grepl('UP', colnames(asv_tab_deblur))]
  asv_tab_deblur = asv_tab_deblur %>% select(-GL140_UP_1_16S_sed, -GL140_UP_2_16S_sed, -GL140_UP_3_16S_sed)
  
  samples_to_keep = colnames(asv_tab_deblur)[colSums(asv_tab_deblur) > min_reads_cutoff]
  asv_tab_deblur = asv_tab_deblur %>% select(any_of(samples_to_keep))
  
  colnames(asv_tab_deblur) = map_chr(colnames(asv_tab_deblur), function(x) paste0(strsplit(x, '_')[[1]][1:2], collapse = '_'))
  return(asv_tab_deblur)}

load_meta_data <- function(asv_tab){
  meta_tab = read.csv('data/nomis-20230320-0855-db.csv')
  meta_tab$sample = map_chr(meta_tab$patch, function(x) paste0(strsplit(x, '_')[[1]][1:2], collapse = '_'))
  meta_tab$rep_n = map_int(meta_tab$sample, function(x) sum(colnames(asv_tab) == x))
  meta_tab = meta_tab %>% select(sample, mountain_range, lat_sn..DD., lon_sp..DD., rep_n) %>% 
    filter(mountain_range != 'Rwenzori Mountains', grepl('_UP', sample)) %>% distinct()
  
  meta_tab = meta_tab[meta_tab$sample %in% colnames(asv_tab),]
  return(meta_tab)}

endemic_analysis <- function(x, meta_tab, n_sample, rarefaction_depth){
  colnames(x) = map_chr(colnames(x), function(y) unique(meta_tab$mountain_range[meta_tab$sample == y]))
  
  # Group by mountain range
  new_x = data.frame(row.names = rownames(x))
  for (mountain_range in unique(colnames(x))){
    samples = sample(meta_tab$sample[meta_tab$mountain_range == mountain_range], n_sample)
    mr_data = x[,colnames(x) == mountain_range]
    new_x = cbind(new_x, mr_data[,sample(ncol(mr_data), n_sample)])}
  colnames(new_x) = map_chr(colnames(new_x), function(x) paste0(strsplit(x, '\\.')[[1]][1]))
  #new_x = as.data.frame(do.call(cbind, by(t(new_x),INDICES=names(new_x),FUN=colSums)))
  new_x = new_x[rowSums(new_x) > 0,]
  
  res = data.frame()
  for (mountain_range in unique(colnames(new_x))){
    endemic_mask = (rowSums(new_x[,colnames(new_x) == mountain_range]) == rowSums(new_x))
    tot_asv_n = sum(rowSums(new_x[,colnames(new_x) == mountain_range]) > 0)
    end_asv_n = sum((rowSums(new_x[,colnames(new_x) == mountain_range]) > 0) & endemic_mask)
    end_asv_median_sample = median(colSums(new_x[endemic_mask, colnames(new_x) == mountain_range] > 0))
    end_asv_median_sample_ab = median(colSums(new_x[endemic_mask, colnames(new_x) == mountain_range]) / colSums(new_x[,colnames(new_x) == mountain_range]))
    res = rbind(res, data.frame(mountain_range=mountain_range, rarefaction_depth=rarefaction_depth, sample_n=n_sample, 
                                tot_asv_n=tot_asv_n, end_asv_n=end_asv_n, 
                                end_asv_prop=end_asv_n/tot_asv_n, end_asv_median_sample_ab=end_asv_median_sample_ab, 
                                end_asv_median_sample_prop=end_asv_median_sample/tot_asv_n,
                                end_asv_median_sample=end_asv_median_sample))}
  return(res)}

endemic_unfiltered_analysis <- function(x, meta_tab, n_sample){
  colnames(x) = map_chr(colnames(x), function(y) unique(meta_tab$mountain_range[meta_tab$sample == y]))
  
  # Group by mountain range
  new_x = data.frame(row.names = rownames(x))
  for (mountain_range in unique(colnames(x))){
    mr_data = x[,colnames(x) == mountain_range]
    new_x = cbind(new_x, mr_data[,sample(ncol(mr_data), n_sample)])}
  colnames(new_x) = map_chr(colnames(new_x), function(x) paste0(strsplit(x, '\\.')[[1]][1]))
  #new_x = as.data.frame(do.call(cbind, by(t(new_x),INDICES=names(new_x),FUN=colSums)))
  new_x = new_x[rowSums(new_x) > 0,]
  
  res = data.frame()
  for (mountain_range in unique(colnames(new_x))){
    endemic_mask = (rowSums(new_x[,colnames(new_x) == mountain_range]) == rowSums(new_x))
    tot_asv_n = sum(rowSums(new_x[,colnames(new_x) == mountain_range]) > 0)
    end_asv_n = sum((rowSums(new_x[,colnames(new_x) == mountain_range]) > 0) & endemic_mask)
    end_asv_median_sample = median(colSums(new_x[endemic_mask, colnames(new_x) == mountain_range] > 0))
    end_asv_median_sample_ab = median(colSums(new_x[endemic_mask, colnames(new_x) == mountain_range]) / colSums(new_x[,colnames(new_x) == mountain_range]))
    res = rbind(res, data.frame(mountain_range=mountain_range, 
                                tot_asv_n=tot_asv_n, end_asv_n=end_asv_n, 
                                end_asv_prop=end_asv_n/tot_asv_n, end_asv_median_sample_ab=end_asv_median_sample_ab, 
                                end_asv_median_sample_prop=end_asv_median_sample/tot_asv_n,
                                end_asv_median_sample=end_asv_median_sample))}
  return(res)}

analyse_endemics <- function(asv_tab_deblur, asv_filtered_tab_deblur, meta_data, n_sample, rd){
  all_res_endemics = data.frame()
  for (i in 1:50){
    set.seed(i)
    # Create summed table, and random replicate table
    tab_filtering_deblur = take_avg(asv_filtered_tab_deblur)
    # Rarefy the tables
    rarefied_tab_filtered_deblur = rarefy_tab(tab_filtering_deblur, rd)
    
    res_filtered_deblur = endemic_analysis(rarefied_tab_filtered_deblur, meta_data, n_sample, rd)
    res_filtered_deblur$Method = 'deblur'
    
    all_res_endemics = rbind(all_res_endemics, res_filtered_deblur)}
  return(all_res_endemics)}

main <- function(){
  asv_tab_deblur = load_asv_data(25000)
  meta_data = load_meta_data(asv_tab_deblur)
upper_bound = endemic_unfiltered_analysis(asv_tab_deblur, meta_data)
  print(upper_bound)
   write.csv(upper_bound, file='endemism_unfiltered.csv', quote=F, row.names=F)
  
  # Custom filtering
  filtered_asvs_deblur = filtering_br(asv_tab_deblur, meta_data)
  asv_filtered_tab_deblur = asv_tab_deblur[rownames(asv_tab_deblur) %in% filtered_asvs_deblur, ]

  
  
  all_res_endemics = data.frame()
  registerDoMC(8)
  
  for (n_samp in c(2,3,4,5,6)){
    print(n_samp)
    res_endemics = foreach(r_depth = c(6580,16580,26580,36580,46580,56580,66580,76580), .combine = 'rbind') %dopar% {analyse_endemics(asv_tab_deblur, asv_filtered_tab_deblur, meta_data, n_samp, r_depth)}
    all_res_endemics = rbind(all_res_endemics, res_endemics)}
  
  
  model_end_n = gam(data=all_res_endemics, family = quasibinomial(),
                        formula = end_asv_prop ~ s(rarefaction_depth, sample_n, by = as.factor(mountain_range), bs='tp', k=3))
      
  n_mr = length(unique(all_res_endemics$mountain_range))
  n_samp = max(all_res_endemics$sample_n)
  n_rar = max(all_res_endemics$rarefaction_depth)
  bounds = data.frame(mountain_range=rep(unique(all_res_endemics$mountain_range), 2), 
                           rarefaction_depth=c(rep(n_rar*2, n_mr), rep(n_rar*5, n_mr)),
                           sample_n=c(rep(n_samp*2, n_mr), rep(n_samp*5, n_mr)),
                           group=c(rep('x2', n_mr), rep('x5', n_mr)))
  bounds$lower_bound = exp(predict(model_end_n, newdata = bounds, se.fit=T)$fit)# - (exp(predict(model_end_n, newdata = bounds, se.fit=T)$se.fit))
  
  
  
  bounds$upper_bound = map_dbl(bounds$mountain_range, function(x) upper_bound$end_asv_prop[upper_bound$mountain_range == x])
}

main()












