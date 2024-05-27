library(phyloseq)
library(phyloseqCompanion)
library(mvabund)

nomis_r <- readRDS("20230825_NOMIS_rarefied.RData")
uganda=c("Uganda")
prune_Uganda <- subset_samples(nomis_r, !site_c %in% uganda)
prune_uganda_metadata <- sample.data.frame(prune_Uganda)

asv_pu <- t(otu_table(prune_Uganda, taxa_are_rows=T))

ab <- mvabund(asv_pu)
asv_nb <- manyglm(ab ~ site_c,
                  data = prune_uganda_metadata, family = 'negative binomial')


nomis_avo <- anova(asv_nb, p.uni= "adjusted", nBoot = 99, show.time=T)
nomis_manyglm_res <- nomis_avo_site$uni.p
write.csv(nomis_manyglm_res, "results_manyglm_nomis.csv")