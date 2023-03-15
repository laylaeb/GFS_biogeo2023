###Shared ASVs entre les paires de régions
## load packages
library(phyloseq)
library(phyloseqCompanion)
library(ggplot2)
library(vegan)
library(dplyr)
library(data.table)
library(DescTools)
library(ComplexUpset)
library("RColorBrewer")
library(speedyseq)
library(ComplexHeatmap)
library(ggforce)
llibrary(tidyverse)

taxglom_sed_family <- tax_glom(physeq=prune_Uganda, taxrank=rank_names(prune_Uganda)[5], NArm=F)
alphadt_NOMIS <- sample.data.frame(prune_Uganda)

Site_c <- c('Nepal.x','Norway.x','Alps.x','Caucasus.x','Ecuador.x','Greenland.x','New_Zealand.x','kyrgyzstan.x','Chile.x','Alaska.x')

alphadt_NOMIS$Site_c <- factor(alphadt_NOMIS$Site_c, levels = Site_c)

ps2.venn <- merge_samples(prune_Uganda, 'Site_c', fun = sum)

venn_obj <- as.data.frame(t(otu_table(ps2.venn)))

## Transform into binary
venn_obj.binary <- sapply(venn_obj, function(x) ifelse(x > 0, 1, 0),
                          USE.NAMES = T)
## Transform relative abundance
venn_obj.ab <- transform_sample_counts(otu_table(venn_obj, taxa_are_rows=T), function(x) x/sum(x))
## Give names to rows
rownames(venn_obj.ab) <- rownames(venn_obj)
rownames(venn_obj.binary) <- rownames(venn_obj)

## Always transform as dataframe
venn_obj_binary_df <- as.data.frame(venn_obj.binary)
venn_obj_ab_df<- as.data.frame(venn_obj.ab)

## Compute Sum of the taxa
venn_obj_ab_df$Ab_Sum <- rowSums(venn_obj_ab_df)

# Merging Table
venn_obj_merge <- merge(venn_obj_binary_df, venn_obj_ab_df, by=0)
venn_obj_merge_df <- as.data.frame(venn_obj_merge)
write.csv(venn_obj_merge_df, "check_shared_percent.csv", sep=",")

sumupset <- as.data.frame(matrix(ncol=3,nrow=0))
colnames(sumupset) <- c("sommec","sommea","n")
for (i in Site_c) {
  for (j in Site_c) {
    if(i==j){next}
    i_y <- gsub(".x",".y",i)
    j_y <- gsub(".x",".y",j)
    i=sym(i)
    j=sym(j)
    i_y=sym(i_y)
    j_y=sym(j_y)
    temp <- venn_obj_merge_df %>% filter(!!i==1 & !!j==1) %>% summarize(sommec=sum(!!i_y),sommea=sum(!!j_y))
    temp2 <- venn_obj_merge_df %>% filter(!!i==1 & !!j==1)  %>% summarize(n=n())
    temp <- as.data.frame(temp)
    name <- paste(i, j, sep="_")
    rownames(temp) <- name
    temp$n <-temp2$n
    sumupset <- rbind(sumupset,temp)
    
  }
}

head(sumupset)

sumupset$shared_mr <- rownames(sumupset)
rownames(sumupset) <- NULL
melt_shared <- melt(sumupset[,c(1,2,4)], id='shared_mr')

ggplot(melt_shared, aes(fill=variable, y=value, x=shared_mr)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip()


##Here sommec and sommea represent the sum of abundance of a mountain range, meaning the contribution of the xx ASVs that are shared in term of relative abundance
##For instance, Nepal share 2157 ASVs with Norway, which account for 41.9% of relative abundance in Nepal. So these 1/4 of ASVs in Nepal contribute to 41.9% of relative abundance
##These are really abundant ones!

## On a pris les resultats de sumupset, on a enlevés les duplicats de lignes (Nepal - Alaska, Alaska - Nepal)
## Ensuite on a mergé les deux colonnes 

percentshared <- read.csv("SharedViolin_2022811.csv")
median(percentshared$percent_pairs)
meltshared <- melt(percentshared)

shared_asv_plot<- ggplot(meltshared, aes(x=variable, y = value*100, fill= variable)) + 
  geom_boxplot(alpha=0.8) + geom_jitter(stat="identity", cex=2, lwd=3)+ 

  theme_bw() +
  ylab("Percent(%)") +
  xlab("")+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),legend.title =element_blank())+
  #scale_fill_manual(labels=c('Relative Abundance','Shared ASVs'),values = wes_palette("Chevalier1"))
  scale_fill_brewer(palette = "Dark2")

shared_asv_plot + coord_flip() + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

##Where proportion ASVs is the proportion of ASVs that are shared between pairs of mountain ranges
##Percent_Pairs represents how much of relative abundance explain those ASVs that are shared 

