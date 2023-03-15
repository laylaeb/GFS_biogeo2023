##Rank abundance Curves, all Glaciers! So we have to take into account the biogeo betadiv --> 127 glaciers
comm<-as(otu_table(prune_Uganda), "matrix")
env<-as(sample.data.frame(prune_Uganda), "matrix")

comm<-as.data.frame(comm)
env<-as.data.frame(env)

### Applying this to the full dataset
comm_melt <- reshape2::melt(as.matrix(comm), na.rm=TRUE)
colnames(comm_melt)<-c("Var1","code_gl","Abundance")
head(comm_melt)
comm_melt<-as.data.frame(comm_melt)
comm_melt_scaling<-comm_melt[comm_melt$Abundance!=0,] ##Removing taxa with 0 abundance
head(comm_melt_scaling)

comm_melt_RA<-comm_melt_scaling %>%
  group_by(code_gl) %>%
  mutate(Scale_RA=Abundance/sum(Abundance)) ##Takes the relative abundance for each sample
head(comm_melt_RA)

##Sanity check
comm_melt_RA %>%
  group_by(code_gl)%>% summarize(somme=sum(Scale_RA))

comm_melt_RANK<-comm_melt_RA %>%
  group_by(code_gl) %>%
  mutate(my_ranks = order(order(Scale_RA,decreasing=T))) #Ranks

## Pour comprendre comment order fonctionne, en fait -- on prendre notre vecteur
## et on regarde à quelle position les chiffres/nombres sont placés et on les replace
## en ordre!
Rank_merge<-merge(comm_melt_RANK, env, by="code_gl")

Rank_merge_mod <- expand.grid(Site_c = unique(Rank_merge$Site_c),
                              Rank = 1:max(Rank_merge$my_ranks))

Rank_merge_mod$Scale_RA = vapply(1:nrow(Rank_merge_mod), 
                                 function(i) mean(Rank_merge$Scale_RA[(Rank_merge$my_ranks == Rank_merge_mod$Rank[i]) & 
                                                                 (Rank_merge$Site_c == Rank_merge_mod$Site_c[i])]),FUN.VALUE = numeric(1))
colors<-c("#2E2A2BFF", "#CF4E9CFF", "#8C57A2FF",
          "#3EBCB6","#82581FFF","#2F509EFF",
          "#E5614CFF","#97A1A7FF","#bee183","#DC9445FF","#EDD03E")
rare_curve <- ggplot() +
  geom_line(data=Rank_merge, aes(x=my_ranks, y=log(Scale_RA), group=code_gl),size=0.2,colour="black",alpha=0.2)+
  geom_line(data=Rank_merge_mod, aes(x=Rank, y=log(Scale_RA), colour=Site_c))+
  scale_color_manual(values= c("Alaska"="#2E2A2BFF","Alps"="#CF4E9CFF","Caucasus"="#8C57A2FF",
                               "Chile"="#3EBCB6","Ecuador"="#82581FFF","Greenland"="#2F509EFF",
                               "Kirghizistan"="#E5614CFF","Nepal"="#97A1A7FF","New_Zealand"="#bee183","Norway"="#DC9445FF"))+
  xlab("Rank")+
  facet_wrap(~Site_c, nrow=3)+
  theme(panel.background = element_rect(fill = 'white', color="grey60"))

## Frequence relative des ASVs
taxa_names(prune_Uganda) <- paste("ASV", 1:ntaxa(prune_Uganda), sep="_")
sort(sample_sums(prune_Uganda))
comm_count <- as.matrix(t(otu_table(prune_Uganda, taxa_are_rows=T)))

spe_table <- data.frame(species=colnames(comm_count),
                        freq.rel=apply(comm_count,2,sum)/sum(comm_count)*100)

spe_table_10 <- subset(spe_table[1:15000,])

setwd(here("Figures"))
pdf("Species_regional_occurence.pdf", width=8, height=15)
ggplot(spe_table, aes(x=reorder(species, freq.rel), y=freq.rel)) +
  geom_bar(position="dodge", stat="identity", fill="steelblue2", colour="grey") +
  xlab("") + ylab("Relative species frequency") +
  coord_flip()+
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()) #remove y axis ticks
dev.off()

library(goeveg)
comm_t <- as.matrix(t(comm))
comm_table_10 <- subset(spe_table[1:5000,])

racu<-racurve(
  comm_t,
  main = "Rank abundance curve",
  nlab = 0,
  ylog = TRUE,
  frequency = FALSE,
  ylim = c(0.0000001, 0.1),
  xlim = NULL
)

###################Character species

##Rank abundance Curves, all Glaciers! So we have to take into account the biogeo betadiv --> 127 glaciers
comm<-as(otu_table(indicator_for_multivariate), "matrix")
env<-as(sample_data(indicator_for_multivariate), "matrix")

comm<-as.data.frame(comm)
env<-as.data.frame(env)

### Applying this to the full dataset
comm_melt <- reshape2::melt(as.matrix(comm), na.rm=TRUE)
colnames(comm_melt)<-c("Var1","code_gl","Abundance")
head(comm_melt)
comm_melt<-as.data.frame(comm_melt)
comm_melt_Scale<-comm_melt[comm_melt$Abundance!=0,] ##Removing taxa with 0 abundance
head(comm_melt_Scale)
comm_melt_RA<-comm_melt_Scale %>%
  group_by(code_gl) %>%
  mutate(Scale_Abundance=Abundance/sum(Abundance)) ##Takes the relative abundance for each sample
head(comm_melt_RA)

##Sanity check
comm_melt_RA %>%
  group_by(code_gl)%>% summarize(somme=sum(Scale_Abundance))

comm_melt_RANK<-comm_melt_RA %>%
  group_by(code_gl) %>%
  mutate(my_ranks = order(order(Scale_Abundance, decreasing=TRUE))) #Ranks

Rank_merge<-merge(comm_melt_RANK, env, by="code_gl")

##Là on prend la moyenne des Scale_abundance par site 
Rank_merge_mod<-Rank_merge %>% group_by(Site_c, Var1) %>% summarize(mean_scale=mean(Scale_Abundance), mean_rank=mean(my_ranks)) 

#Getting the plot
ggplot(data=Rank_merge_mod, aes(x=mean_rank, y=mean_scale, color=Site_c)) +
  #geom_point(alpha=0.1)+
  geom_smooth(method="loess",se=F)+
  #geom_line(size=0.5)+
  xlab("Rank")+
  ylab("Relative Abundance")+
  #scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')+
  ggtitle("All glaciers")+
  facet_wrap(~Site_c, nrow=3)+
  theme_bw()

#Getting the plot
ggplot(data=Rank_merge, aes(x=my_ranks, y=Scale_Abundance, group=code_gl, color=Site_c)) +
  #geom_point(alpha=0.1)+
  #geom_smooth(method="loess",se=F)+
  geom_line(size=0.5)+
  xlab("Rank")+
  ylab("Relative Abundance")+
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')+
  ggtitle("All glaciers")+
  facet_wrap(~Site_c, nrow=3)+
  #ylim(log10(0), log10(0.5))+
  theme_bw()


## Frequence relative des ASVs
taxa_names(indicator_for_multivariate) <- paste("ASV", 1:ntaxa(indicator_for_multivariate), sep="_")
sort(sample_sums(indicator_for_multivariate))
comm_count <- as.matrix(t(otu_table(indicator_for_multivariate, taxa_are_rows=T)))
spe_table <- data.frame(species=colnames(comm_count),
                        freq.rel=apply(comm_count,2,sum)/sum(comm_count)*100)

spe_table_10 <- subset(spe_table[1:150,])

setwd(here("Figures"))
pdf("Species_regional_occurence.pdf", width=8, height=15)
ggplot(spe_table_10, aes(x=reorder(species, freq.rel), y=freq.rel)) +
  geom_bar(position="dodge", stat="identity", fill="steelblue2", colour="grey") +
  xlab("") + ylab("Relative species frequency") +
 # coord_flip()+
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()) #remove y axis ticks
dev.off()

library(goeveg)
comm_t <- as.matrix(t(comm))
comm_table_10 <- subset(spe_table[1:5000,])

racu<-racurve(
  comm_t,
  main = "Rank abundance curve",
  nlab = 0,
  ylog = TRUE,
  frequency = FALSE,
  ylim = c(0.0000001, 0.1),
  xlim = NULL
)

### Autre essai
env$Site_c <- as.factor(env$Site_c)
essai2 <- rankabundance((comm),  digits=1, t=qt(0.975, df=n-1))
