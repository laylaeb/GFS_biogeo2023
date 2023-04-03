setwd(choose.dir())
ASVs<-read.csv("20221026NOMIS_full_asv_sedimentup_filtered_final.csv", head=T, row.names = 1)
ASVs.unfiltered<-read.csv("20230103_NOMIS_unfiltered_sedimentup_table.csv", sep="",row.names = 1)
tax<-read.csv("20221026NOMIS_full_tax_sedimentup_filtered_final.csv", head=T, row.names=1)
ASVs.endemic<-read.csv("endemic_asv_table.csv")
samp<-read.csv("samples148.csv")



# taxonomic composition ---------------------------------------------------

library(metacoder)

tax<-read.csv("20221026NOMIS_taxtree.csv",head=T, row.names=1)

obj_v <- parse_tax_data(tax, class_cols="lineage", class_sep = ";")
obj_v$data$tax_abund <- calc_taxon_abund(obj_v, "tax_data")
obj_v$data$tax_occ <- calc_n_samples(obj_v, "tax_abund")

heat_tree(obj_v, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = obj_v$data$tax_occ$n_samples,
          node_label_size=n_samples,
          node_label_size_range = c(0.005, 0.01),
          node_color_axis_label = "occupancy",
          node_size_axis_label="",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford")



