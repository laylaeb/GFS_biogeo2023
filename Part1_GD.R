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

uganda=c("Uganda")
prune_Uganda <- subset_samples(NOMIS_FR, !Site_c %in% uganda)

site_c_g = c("Greenland")
greenland <- subset_samples(prune_Uganda, Site_c %in% site_c_g)
glid_sample_greenland <- merge_samples(greenland, "code_gl")

site_c_no = c("Norway")
norway<- subset_samples(prune_Uganda, Site_c %in% site_c_no)
glid_sample_norway <- merge_samples(norway, "code_gl")

site_c_al = c("Alps")
alps <- subset_samples(prune_Uganda, Site_c %in% site_c_al)
glid_sample_alps <- merge_samples(alps, "code_gl")

site_c_ec = c("Ecuador")
ecuador <- subset_samples(prune_Uganda, Site_c %in% site_c_ec)
glid_sample_ecuador<- merge_samples(ecuador, "code_gl")

site_c_ne = c("Nepal")
nepal <- subset_samples(prune_Uganda, Site_c %in% site_c_ne)
glid_sample_nepal<- merge_samples(nepal, "code_gl")

site_c_ca = c("Caucasus")
caucasus <- subset_samples(prune_Uganda, Site_c %in% site_c_ca)
glid_sample_caucasus<- merge_samples(caucasus, "code_gl")

site_c_nz= c("New_Zealand")
nz <- subset_samples(prune_Uganda, Site_c %in% site_c_nz)
glid_sample_nz<- merge_samples(nz, "code_gl")

site_c_k = c("Kirghizistan")
kh <- subset_samples(prune_Uganda, Site_c %in% site_c_k)
glid_sample_kh<- merge_samples(kh, "code_gl")

site_c_chile = c("Chile")
chile <- subset_samples(prune_Uganda, Site_c %in% site_c_chile)
glid_sample_chile<- merge_samples(chile, "code_gl")

site_c_alaska = c("Alaska")
alaska <- subset_samples(prune_Uganda, Site_c %in% site_c_alaska)
glid_sample_alaska<- merge_samples(alaska, "code_gl")

#Alps
asv_spcom_alps <- otu_table(glid_sample_alps, taxa_are_rows = T)
asv_spcom_alps_t <- t(asv_spcom_alps)

#Caucasus
asv_spcom_caucasus<- otu_table(glid_sample_caucasus, taxa_are_rows = T)
asv_spcom_caucasus_t <- t(asv_spcom_caucasus)

#Ecuador
asv_spcom_ecuador<- otu_table(glid_sample_ecuador, taxa_are_rows = T)
asv_spcom_ecuador_t <- t(asv_spcom_ecuador)

#Greenland
asv_spcom_greenland<- otu_table(glid_sample_greenland, taxa_are_rows = T)
asv_spcom_greenland_t<- t(asv_spcom_greenland)

#NZ
asv_spcom_nz<- otu_table(glid_sample_nz, taxa_are_rows = T)
asv_spcom_nz_t <- t(asv_spcom_nz)

#Nepal
asv_spcom_kh<- otu_table(glid_sample_nepal, taxa_are_rows = T)
asv_spcom_kh_t <- t(asv_spcom_kh)

#Kirghizistan
asv_spcom_nepal<- otu_table(glid_sample_kh, taxa_are_rows = T)
asv_spcom_nepal_t <- t(asv_spcom_nepal)

#Norway
asv_spcom_norway<- otu_table(glid_sample_norway, taxa_are_rows = T)
asv_spcom_norway_t <- t(asv_spcom_norway)

#chile
asv_spcom_chile<- otu_table(glid_sample_chile, taxa_are_rows = T)
asv_spcom_chile_t <- t(asv_spcom_chile)

#Alaska
asv_spcom_alaska<- otu_table(glid_sample_alaska, taxa_are_rows = T)
asv_spcom_alaska_t <- t(asv_spcom_alaska)

##Prepare frequency table for iNEXT
##Greenland
asv_region_greenland_df <- asv_spcom_greenland_t
sumrow_greenland <- unname(rowSums(asv_region_greenland_df>0))
sort_greenland <- sort(sumrow_greenland, decreasing=T)
vec_greenland <- sort_greenland[sort_greenland >0]

##Norway
asv_region_norway_df <- asv_spcom_norway_t
sumrow_norway <- unname(rowSums(asv_region_norway_df>0))
sort_norway<- sort(sumrow_norway, decreasing=T)
vec_norway <- sort_norway[sort_norway >0]

##Ecuador
asv_region_ecuador_df <- asv_spcom_ecuador_t
sumrow_ecuador <- unname(rowSums(asv_region_ecuador_df>0))
sort_ecuador<- sort(sumrow_ecuador, decreasing=T)
vec_ecuador <- sort_ecuador[sort_ecuador >0]

##Alps
asv_region_alps_df <- asv_spcom_alps_t
sumrow_alps <- unname(rowSums(asv_region_alps_df>0))
sort_alps<- sort(sumrow_alps, decreasing=T)
vec_alps <- sort_alps[sort_alps >0]

##Nepal
asv_region_nepal_df <- asv_spcom_nepal_t
sumrow_nepal <- unname(rowSums(asv_region_nepal_df>0))
sort_nepal<- sort(sumrow_nepal, decreasing=T)
vec_nepal <- sort_nepal[sort_nepal >0]

##KH
asv_region_kh_df <- asv_spcom_kh_t
sumrow_kh <- unname(rowSums(asv_region_kh_df>0))
sort_kh<- sort(sumrow_kh, decreasing=T)
vec_kh <- sort_kh[sort_kh >0]

##Caucasus
asv_region_caucasus_df <- asv_spcom_caucasus_t
sumrow_caucasus <- unname(rowSums(asv_region_caucasus_df>0))
sort_caucasus<- sort(sumrow_caucasus, decreasing=T)
vec_caucasus <- sort_caucasus[sort_caucasus >0]

##NZ
asv_region_nz_df <- asv_spcom_nz_t
sumrow_nz <- unname(rowSums(asv_region_nz_df>0))
sort_nz<- sort(sumrow_nz, decreasing=T)
vec_nz <- sort_nz[sort_nz >0]

##chile
asv_region_chile_df <- asv_spcom_chile_t
sumrow_chile <- unname(rowSums(asv_region_chile_df>0))
sort_chile<- sort(sumrow_chile, decreasing=T)
vec_chile <- sort_chile[sort_chile >0]

##alaska
asv_region_alaska_df <- asv_spcom_alaska_t
sumrow_alaska <- unname(rowSums(asv_region_alaska_df>0))
sort_alaska<- sort(sumrow_alaska, decreasing=T)
vec_alaska <- sort_alaska[sort_alaska >0]

library(iNEXT)
list_exped_all <- list(alps=c(ncol(asv_region_alps_df),vec_alps),nz=c(ncol(asv_region_nz_df),vec_nz),
                       caucasus=c(ncol(asv_region_caucasus_df),vec_caucasus),
                       kh=c(ncol(asv_region_kh_df),vec_kh),
                       greenland=c(ncol(asv_region_greenland_df),vec_greenland),
                       norway=c(ncol(asv_region_norway_df),vec_norway), 
                       ecuador=c(ncol(asv_region_ecuador_df),vec_ecuador),
                       nepal=c(ncol(asv_region_nepal_df),vec_nepal),
                       alaska=c(ncol(asv_region_alaska_df),vec_alaska),
                       chile=c(ncol(asv_region_alaska_df),vec_chile))
                       

t <- seq(1, 30, by=5)
out_all_exped <- iNEXT(list_exped_all, q=0, datatype="incidence_freq", size=t)
