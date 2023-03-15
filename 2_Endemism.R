##Some metrics regarding the Number of ASVs found across the dataset

NOMIS_FR <- readRDS("NOMIS_FR_20221129.RData")
uganda=c("Uganda")
prune_Uganda <- subset_samples(NOMIS_FR, !Site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))

metadata_nomis <- sample.data.frame(prune_Uganda)
asv_df <- as.data.frame(t(otu_table(prune_Uganda, taxa_are_rows=T)))

## Here we are investigating the prevalence of ASV across the dataset
howmanyasv<- as.data.frame(colSums(asv_df != 0))
colnames(howmanyasv) <- c("Count_nb")
howmanyasv$ASV <- rownames(howmanyasv)
rownames(howmanyasv) <- NULL

## Endemic ASVs
avsdfmelt <- melt(as.matrix(asv_df))
## Keep only values that are > 0
avsdfmelt <- avsdfmelt[avsdfmelt$value >0,]

##### Endemism
## Liste des ASVs qui ont une prévalence de 1 : 18492 soit 40.8% -- comment sont-ils répartis par mountain ranges?
equalone <- howmanyasv %>% filter(Count_nb == 1)
## On prend notre table d'ASV et on filtre pour avoir le nom des échantillons ## value c'est le nb de count, et count_nb c'est le nb d'échantillons
merge_equalone <- merge(equalone, avsdfmelt, by.x="ASV",by.y="Var2")
merge_equalone$MR <- vapply((merge_equalone$Var1), function(x) metadata_nomis$Site_c[metadata_nomis$code_gl == x], FUN.VALUE = character(1))
colnames(merge_equalone) <- c("ASV","prev","glname","nb_count","MR")
##Mtn on va obtenir le nombre d'ASVs qui sont endemiques par région -- prenant en compte une prévalence de 1 (observé 1x)
essai_final <- merge_equalone %>% group_by(MR) %>% summarize(prev=sum(prev)) ##18492
essai_plot <- merge_equalone %>% group_by(MR) 
essai_plotv2 <- essai_plot[c("MR","ASV","Count_nb")]
colnames(essai_plotv2) <- c("MR","ASV","prev")
essai_plotv2$Color <- "C_uniquetoone"

# MR              prev
# <chr>          <dbl>
# 1 Alaska        1546 -> 1546/18492= 8.4
# 2 Alps          2261 -> 2261/18492 =12.2
# 3 Caucasus      1650 -> 1650/18492 = 8.9
# 4 Chile         2055 -> 2055/18492 = 11.1
# 5 Ecuador       2012 -> 2012/18492 = 10.9
# 6 Greenland      452 -> 452/18492 = 2.4
# 7 kyrgyzstan  1546 -> 1546/18492 = 8.4
# 8 Nepal         2617 -> 2617/18492 = 14.2
# 9 New_Zealand   2555 -> 2555/18492 = 13.8
# 10 Norway       1798 -> 1798/18492 = 9.7


## Barplot 
df_equalone <- data.frame(Mountain_Range=c("Alaska","European Alps","Caucasus Mountains","Chilean Andes","Ecuadorian Andes","Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps",
                                             "Scandinavian Mountains"), Percentage=c(8.4, 12.2, 8.9, 11.1, 10.9, 2.4,8.4,14.2,13.8,9.7))
## Reorder
df_equalone$Mountain_Range<- factor(df_equalone$Mountain_Range, levels = df_equalone$Mountain_Range[order(df_equalone$Percentage)])

p<-ggplot(data=df_equalone, aes(x=Mountain_Range, y=Percentage)) +geom_bar(stat="identity")
p + coord_flip() + theme_minimal()


## Filter the table of prevalence 
twoandnine<- howmanyasv %>% filter(Count_nb > 1 & Count_nb < 10) ## 21386 ASVs --47.1% des ASVs
## On merge twoandnine avec la table d'ASV original pour obtenir le nom des échantillons puis pour chaque échantillon on lui attribue un nom de MR!
merge_twoandnine <- merge(twoandnine, avsdfmelt, by.x="ASV",by.y="Var2")
merge_twoandnine$MR <- vapply((merge_twoandnine$Var1), function(x) metadata_nomis$Site_c[metadata_nomis$code_gl == x], FUN.VALUE = character(1))
## Maintenant on va compter le nombre de fois que l'ASV est présent dans la region, donc on va summer par GL code!(et on somme la value)
twoandnine_end<- merge_twoandnine %>% group_by(MR,ASV) %>% summarize(prev=n())
twoandnine_end$Color <- "B_twoandnine"
## Maintenant on veut qu'il soit unique dans la mountain range! On  mutate pour pouvoir compter le nombre de fois ou l'ASV est observé, puis on filtre seulement les ASVs qui sont observés 1X dans la MR
## Après cela on va grouper par MR et on va faire la sum des N pour savoir combien il y a d'ASVs dans chaque MR
endemism_twoandnine <- twoandnine_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n()) ## 7779 ASVs out of 21386 soit 36.4 %

# A tibble: 10 × 2
# MR           number
# <chr>         <int>
# 1 Alaska           285 -- 285/7779 =3.66
# 2 Alps            908 -- 908/7779= 11.67
# 3 Caucasus        557 -- 557/7779 = 7.16
# 4 Chile           492 -- 492/7779 = 6.32
# 5 Ecuador         830 -- 830/7779 = 10.7
# 6 Greenland        65 -- 65/7779 = 0.84
# 7 kyrgyzstan    879 -- 879/7779 = 11.29
# 8 Nepal          1038 -- 1038/7779 = 13.34
# 9 New_Zealand    2346 -- 2346/7779 = 30.16
# 10 Norway         379 -- 379/7779 = 4.87

endemism_twoandnine_plot <- twoandnine_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)

## Barplot 
df_twoandnine <- data.frame(Mountain_Range=c("Alaska","European Alps","Caucasus Mountains","Chilean Andes","Ecuadorian Andes","Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps",
                                             "Scandinavian Mountains"), Percentage=c(3.66, 11.67, 7.16, 6.32, 10.7, 0.84,11.29,13.34,30.16,4.87))
## Reorder
df_twoandnine$Mountain_Range<- factor(df_twoandnine$Mountain_Range, levels = df_twoandnine$Mountain_Range[order(df_twoandnine$Percentage)])

p<-ggplot(data=df_twoandnine, aes(x=Mountain_Range, y=Percentage)) +geom_bar(stat="identity")
p + coord_flip() + theme_minimal()
                            
## On peut faire pareil avec la dernière tranche 
tenandmore<- howmanyasv %>% filter(Count_nb >= 10) ## 5492 ASVs
## On merge twoandnine avec la table d'ASV original pour obtenir le nom des échantillons puis pour chaque échantillon on lui attribue un nom de MR!
merge_tenandmore <- merge(tenandmore, avsdfmelt, by.x="ASV",by.y="Var2")
merge_tenandmore$MR <- vapply((merge_tenandmore$Var1), function(x) metadata_nomis$Site_c[metadata_nomis$code_gl == x], FUN.VALUE = character(1))
tenandmore_end<- merge_tenandmore %>% group_by(MR,ASV) %>% summarize(prev=sum(value>0))

tenandmore_end$Color <-"A_tenandemore"
## Maintenant on veut qu'il soit unique dans la mountain range! On  mutate pour pouvoir compte le nombre de fois ou l'ASV est observé, puis on filtre seulement les ASVs qui sont observés 1X dans la MR
## Après cela on va grouper par MR et on va faire la sum des N pour savoir combien il y a d'ASVs dans chaque MR
endemism_tenandmore <- tenandmore_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n()) ## 161 ASVs out of 5492 2.9%

endemism_tenandmore_plot <- tenandmore_end %>% group_by(ASV) %>% mutate(n=n()) %>% filter(n==1)

# # A tibble: 5 × 2
# MR           number
# <chr>         <int>
#1 Alps            6 -- 6/161 = 3.73
# 2 Caucasus          6 -- 6/161= 3.73
# 3 kyrgyzstan     37 -- 37 / 161 =22.98
# 4 Nepal             2 -- 2/161=1.24
# 5 New_Zealand     110 -- 110/161 = 68.32

df_tenandmore<- data.frame(Mountain_Range=c("European Alps","Caucasus Mountains","Pamir & Tien Shan","Himalayas","Southern Alps"), 
                           Percentage=c(3.73, 3.73, 22.98, 1.24, 68.32))
## Reorder
df_tenandmore$Mountain_Range<- factor(df_tenandmore$Mountain_Range, levels = df_tenandmore$Mountain_Range[order(df_tenandmore$Percentage)])

p<-ggplot(data=df_tenandmore, aes(x=Mountain_Range, y=Percentage)) +geom_bar(stat="identity")
p + coord_flip() + theme_minimal()

## Control
control_asv <- avsdfmelt
control_asv$MR <- vapply((control_asv$Var1), function(x) metadata_nomis$Site_c[metadata_nomis$code_gl == x], FUN.VALUE = character(1))

controlasv_end<- control_asv %>% group_by(MR,Var2) %>% summarize(prev=n()) %>%
  ungroup()%>% group_by(Var2) %>% mutate(n=n())%>% filter(n==1)%>%
  ungroup()%>% group_by(MR) %>%summarize(number=n())

##26432
# 
# # A tibble: 10 × 2
# MR           number
# <chr>         <int>
# 1 Alaska         1831 -> 1831/45370 =4.0        1831/7072=25.9%
# 2 Alps           3175 -> 3175/45370 =7          3175/14130=22.5%
# 3 Caucasus       2213 -> 2213/45370 =4.9        2213/12069=18.3%
# 4 Chile          2547 -> 2547/45370 =5.6        2547/7583=33.6%
# 5 Ecuador        2842 -> 2842/45370 =6.3        2842/7121=39.9%
# 6 Greenland       517 -> 517/45370 = 1.1        517/4033=12.8%
# 7 kyrgyzstan   2462 -> 2462/45370= 5.4        2462/8515=28.9%
# 8 Nepal          3657 -> 3657/45370= 8.1        3657/12332=29.7%
# 9 New_Zealand    5011 -> 5022/45370= 11.1       5022/10802=46.5%
# 10 Norway        2177 -> 2177/45370= 4.8        2177/7764=28%

## to get to know how many ASVs are present per mountain range (Observed)
controlasv_total<- control_asv %>% group_by(MR,Var2) %>% summarize(sumi=sum(n()))%>%
  ungroup()%>% group_by(MR)%>%summarize(sumii=sum(n()))

# MR           sumii
# <chr>        <int>
#   1 Alaska     7072
# 2 Alps         14130
# 3 Caucasus     12069
# 4 Chile         7583
# 5 Ecuador       7121
# 6 Greenland     4033
# 7 kyrgyzstan  8515
# 8 Nepal        12332
# 9 New_Zealand  10802
# 10 Norway       7764

##Now we rbind the 3 different dataframes to plot the ASVs that are endemic to the different mountainranges
essai_plotv2$n <- 1
df_full <- rbind(essai_plotv2, endemism_tenandmore_plot, endemism_twoandnine_plot)
niveaux <- c("New_Zealand","Nepal","Alps","Ecuador","Chile","kyrgyzstan","Caucasus","Norway","Alaska","Greenland")
df_full$MR<- factor(df_full$MR, levels = niveaux)
df_full$MR <- fct_rev(df_full$MR)
write.csv(df_full, "endemic_asv_table.csv")
df_full <- read.csv("endemic_asv_table.csv", sep=",", header=T)
## Plot Stackplot
ggplot(df_full, aes(fill=Color, y=MR)) + 
  geom_bar(position="stack", stat="count")+theme_minimal()

## We would like create a stackedplot with the taxonomy of the 3 different Fractions 
## For all the taxa but also for the endemic only! =)
## Lets do for all taxa from blue/salmon/beige. 

merge_asv_endemic<- merge(t(asv_df), df_full, by.x="row.names",by.y="ASV")
##remove the columns that we dont need anymore
merge_asv_endemic<- as.data.frame(merge_asv_endemic[,-c(149:152)])
##create phyloseq object
row.names(merge_asv_endemic)<-merge_asv_endemic$Row.names
merge_asv_endemic$Row.names <- NULL

endemic_table <- otu_table(merge_asv_endemic, taxa_are_rows=T)
write.csv(endemic_table,"endemic_table_20221412.csv")
merge_endemic_phylo <- merge_phyloseq(endemic_table, tax_NOMIS, metadata_glaciers)
saveRDS(merge_endemic_phylo, "merge_endemic_phylo.RDS")
             
endemic_taxglom <- tax_glom(merge_endemic_phylo, taxrank=rank_names(merge_endemic_phylo)[5], NArm=F)
transf_endemic = transform_sample_counts(endemic_taxglom, function(x) x / sum(x))
sample_merge_region = merge_samples(transf_endemic, "Site_c")
region_endemic = transform_sample_counts(sample_merge_region, function(x) x / sum(x))

TopASV_f <- names(sort(taxa_sums(region_endemic), TRUE)[1:19])
top15_NOMIS_f <- prune_species(TopASV_f, region_endemic)
top15_NOMIS_f <- prune_taxa(taxa_sums(top15_NOMIS_f)>0, top15_NOMIS_f)
top_family<-as.data.frame(tax_table(top15_NOMIS_f))

##Turn all ASVs into Family counts
endemic_df <- psmelt(region_endemic) # create dataframe from phyloseq object
endemic_df$Family<- as.character(endemic_df$Family) #convert to character

## On va choisir de mettre en Other uniquement les ASVs qui ne font pas partie des Familles les plus abondantes
endemic_df$Family[!(endemic_df$Family %in% top_family$Family)] <- "Other"
endemic_df$Family[(endemic_df$Family == "")] <- "Other"
endemic_df$Family[(endemic_df$Family == "g__uncultured")] <- "Other"

##Generating n colors via rbrewer
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

barplot_biogeo <- ggplot(data=endemic_df, aes(x=Sample, y=Abundance, fill=Family))
barplot_biogeo + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A","#66A61E",
                               "#E6AB02","#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                               "#FFFF99")) +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))


## % d'abondance relative des ASVS.
#1) prendre l'abondance relative de toute la table d'ASVs
#2) merger la table avec la liste des endemiques
#3) faire un object phyloseq...et filtrer

prune_uganda_abondance= transform_sample_counts(prune_Uganda, function(x) x / sum(x))
prune_uganda_asv <- otu_table(prune_uganda_abondance, taxa_are_rows=T)
prune_uganda_asv <- prune_uganda_asv[rowSums(prune_uganda_asv[])>0,]
merge_endemic_abondance <- merge(df_full, prune_uganda_asv, by.x="ASV", by.y="row.names")

##remove the columns that we dont need anymore
merge_endemic_abondance<- as.data.frame(merge_endemic_abondance[,-c(2:5)])
##create phyloseq object
row.names(merge_endemic_abondance)<-merge_endemic_abondance$ASV
merge_endemic_abondance$ASV <- NULL
endemic_table_abondance <- otu_table(merge_endemic_abondance, taxa_are_rows=T)
merge_endemic_abondance_phylo <- merge_phyloseq(endemic_table_abondance, tax_NOMIS,metadata_glaciers)

end_ab_phylo_table <- (as.matrix(otu_table(merge_endemic_abondance_phylo, taxa_are_rows=T)))

melt_asv <- melt(end_ab_phylo_table)
merge_asv_data <- merge(as.data.frame(melt_asv), as.matrix(metadata_glaciers), by.x="Var2",by.y="code_gl")

##Here we would need to divise by the number of glaciers per mountain ranges. 
sum_mr <- merge_asv_data %>% group_by(Site_c)%>% summarize(summrr=sum(value), n=n_distinct(Var2))%>% summarize(ar_mr=summrr/n, Site_c)

# 1 0.0961 Alaska      
# 2 0.0437 Alps        
# 3 0.0666 Caucasus    
# 4 0.226  Chile       
# 5 0.253  Ecuador     
# 6 0.0261 Greenland   
# 7 0.146  kyrgyzstan
# 8 0.111  Nepal       
# 9 0.301  New_Zealand 
# 10 0.0577 Norway  
# 
# mean(sum_mr$ar_mr)
# [1] 0.1327927
# > sd(sum_mr$ar_mr)
# [1] 0.09595693


## The most abundant phyla for the endemic ASVs

### First we tax_glom the dataset
#merge_endemic_phylo <- merge_phyloseq(endemic_table, tax_NOMIS, metadata_nomis)

endemic_taxglom_phyla <- tax_glom(merge_endemic_phylo, taxrank=rank_names(merge_endemic_phylo)[2], NArm=F)
tax_table_end_phyla <- tax_table(endemic_taxglom_phyla)
transf_endemic_phyla = transform_sample_counts(endemic_taxglom_phyla, function(x) x / sum(x))
endemic_region = merge_samples(transf_endemic_phyla, "Site_c")
trans_ra_end_phy = transform_sample_counts(endemic_region, function(x) x / sum(x))

asv_endemic_phy <- otu_table(trans_ra_end_phy, taxa_are_rows=T)
melt_endemic_phy <- psmelt(asv_endemic_phy)
merge_taxo_endemic <- merge(melt_endemic_phy, tax_table_end_phyla, by.x="OTU", by.y="row.names")

##Now we would like to know what are the most abundant phyla across the dataset.
sumtot_endo_phylum <- merge_taxo_endemic %>% group_by(Phylum) %>% summarize(sum=sum(Abundance/10))%>%
  filter(!(Phylum %in% c(""," g__uncultured"))) %>%
  slice_max(n=15, order_by=sum)

# A tibble: 15 × 2
# Phylum                      sum
# <chr>                     <dbl>
#   1 " p__Proteobacteria"  0.476  
# 2 " p__Bacteroidota"      0.154  
# 3 " p__Verrucomicrobiota" 0.0705 
# 4 " p__Patescibacteria"   0.0609 
# 5 " p__Planctomycetota"   0.0418 
# 6 " p__Acidobacteriota"   0.0290 
# 7 " p__Actinobacteriota"  0.0236 
# 8 " p__Bdellovibrionota"  0.0207 
# 9 " p__Chloroflexi"       0.0185 
# 10 " p__Myxococcota"       0.0177 
# 11 " p__Gemmatimonadota"   0.0176 
# 12 " p__Cyanobacteria"     0.0141 
# 13 " p__Nitrospirota"      0.0125 
# 14 " p__Armatimonadota"    0.00760
# 15 " p__Desulfobacterota"  0.00482


## Most abundant genera

endemic_taxglom_genera <- tax_glom(merge_endemic_phylo, taxrank=rank_names(merge_endemic_phylo)[6], NArm=F)
tax_table_end_genera <- tax_table(endemic_taxglom_genera)
transf_endemic_genera = transform_sample_counts(endemic_taxglom_genera, function(x) x / sum(x))
endemic_region_genera = merge_samples(transf_endemic_genera, "Site_c")
trans_ra_end_gen = transform_sample_counts(endemic_region_genera, function(x) x / sum(x))

asv_endemic_gen <- otu_table(trans_ra_end_gen, taxa_are_rows=T)
melt_endemic_gen <- psmelt(asv_endemic_gen)
merge_taxo_endemic_gen <- merge(melt_endemic_gen, tax_table_end_genera, by.x="OTU", by.y="row.names")

##Now we would like to know what are the most abundant genera across the dataset.
sumtot_endo_genera <- merge_taxo_endemic_gen %>% group_by(Genus) %>% summarize(sum=sum(Abundance/10))%>%
  filter(!(Genus %in% c(""," g__uncultured"))) %>%
  slice_max(n=15, order_by=sum)

# A tibble: 15 × 2
# Genus                              sum
# <chr>                            <dbl>
# 1 " g__Rhodoferax"                0.0509

# 2 " g__Methylotenera"             0.0467
# 3 " g__Flavobacterium"            0.0390
# 4 " g__Polaromonas"               0.0268
# 5 " g__Rhizobacter"               0.0259
# 6 " g__Candidatus_Nitrotoga"      0.0202
# 7 " g__Candidatus_Nomurabacteria" 0.0198
# 8 " g__Sideroxydans"              0.0170
# 9 " g__Thiobacillus"              0.0165
# 10 " g__Chthoniobacter"            0.0161
# 11 " g__Ferruginibacter"           0.0143
# 12 " g__Ellin6067"                 0.0137
# 13 " g__Hymenobacter"              0.0127
# 14 " g__Gallionella"               0.0120
# 15 " g__0319-6G20"                 0.0107


## How many phylum, classes, orders, families and genera do we have??
phylumfac = factor(tax_table(prune_Uganda)[, "Phylum"])
classfac = factor(tax_table(prune_Uganda)[, "Class"])
orderfac = factor(tax_table(prune_Uganda)[, "Order"])
familyfac = factor(tax_table(prune_Uganda)[, "Family"])
genusfac = factor(tax_table(prune_Uganda)[, "Genus"])
#47 phyla
#143 Class
#Order 373
#Family 564
#Genus 1004

