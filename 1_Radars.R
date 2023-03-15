###Computing radarplots from environmental dataset
require(ggplot2)
require(reshape)
require(scales)
library(scales)
library(ggthemes)
library(ggpubr)

uganda=c("Uganda")
prune_Uganda <- subset_samples(NOMIS_FR, !Site_c %in% uganda)
prune_Uganda_df <- as.matrix(otu_table(prune_Uganda, taxa_are_rows=T))
metadata_nomis_pca <- sample.data.frame(prune_Uganda)

##Variables that we will keep: temperature, pH, conductivity, turbidity, DOC, chla, gl_sa and gl_cov, DIN and SRP
metadata_nomis_pca$DIN <- as.numeric(metadata_nomis_pca$DIN)
metadata_nomis_pca$srp <- as.numeric(metadata_nomis_pca$srp)
metadata_nomis_pca$turb <- as.numeric(metadata_nomis_pca$turb)
metadata_nomis_pca$doc <- as.numeric(metadata_nomis_pca$doc)
metadata_nomis_pca$Chla <- as.numeric(metadata_nomis_pca$Chla)
metadata_nomis_pca$Clays <- as.numeric(metadata_nomis_pca$Clays)
metadata_nomis_pca$Feldspar <- as.numeric(metadata_nomis_pca$Feldspar)
metadata_nomis_pca$Quartz <- as.numeric(metadata_nomis_pca$Quartz)
metadata_nomis_pca$Calcite <- as.numeric(metadata_nomis_pca$Calcite)
metadata_nomis_pca$lon_sp <- abs(metadata_nomis_pca$lon_sp)
metadata_nomis_pca$lat_sp <- abs(metadata_nomis_pca$lat_sp)
metadata_nomis_pca_subset <- metadata_nomis_pca[c('Site_c','water_temp','pH','conductivity','turb','gl_sa','gl_cov','srp','DIN','doc','Chla','Clays','Feldspar','Quartz','Calcite')]
## Remove nas 
df_na <-(na.omit(metadata_nomis_pca_subset))
## Check data distribution
## Histogram 
library(Hmisc)
hist.data.frame(df_na)
log_const <- function(x){return(log(x + (min(x[which(x > 0)])/2)))}
nomis_wt <- log_const(df_na$water_temp)
hist(nomis_wt)
nomis_ph <- log_const(df_na$pH)
hist(nomis_ph)
nomis_cond <- log_const(df_na$conductivity)
hist(nomis_cond)
nomis_turb <- log_const(df_na$turb)
hist(nomis_turb)
nomis_doc <- log_const(df_na$doc)
hist(nomis_doc)
nomis_DIN <- log_const(df_na$DIN)
hist(nomis_DIN)
nomis_srp <- log_const(df_na$srp)
hist(nomis_srp)
nomis_chl <- log_const(df_na$Chla)
hist(nomis_chl)
nomis_glsa <- log_const(df_na$gl_sa)
hist(nomis_glsa)
nomis_Quartz <- log_const(df_na$Quartz)
hist(nomis_Quartz)
nomis_Calcite<- log_const(df_na$Calcite)
hist(nomis_Calcite)
nomis_feldspar <- log_const(df_na$Feldspar)
hist(nomis_feldspar) 
nomis_clays<-asin(df_na$Clays)
hist(nomis_clays)
nomis_glcov <- df_na$gl_cov

df_na_transformed <- as.data.frame(as.matrix(t(rbind(df_na$Site_c,nomis_wt, nomis_glcov, nomis_ph, nomis_cond, nomis_turb, nomis_doc, nomis_DIN, nomis_srp, nomis_chl, nomis_glsa, nomis_DIN, nomis_Quartz, nomis_feldspar, nomis_clays))))
df_num <- as.data.frame(sapply(df_na_transformed[,2:15], as.numeric)) 
df_num_transformed <- cbind(df_num, Site_c=df_na_transformed$V1)
metadata_meansd <- as.data.frame(df_num_transformed %>% group_by(Site_c) %>% summarize_if(is.numeric, list("mean"=mean,"sd"=sd,"min"=min,"max"=max)))

metadata_temp <- as.data.frame(metadata_meansd[, grepl("wt",colnames(metadata_meansd))])
rownames(metadata_temp) <- c("Alaska","Alps","Caucasus","Chile","Ecuador","Greenland","kyrgyzstan","Nepal","New_Zealand","Norway")
metadata_temp$Site_c <- rownames(metadata_temp)
rownames(metadata_temp) <- NULL
head(metadata_temp)

# nomis_wt_mean nomis_wt_sd nomis_wt_min nomis_wt_max    Site_c
# 1    -0.4106212   0.9944335    -1.622017     1.162369    Alaska
# 2    -0.6424682   0.9743456    -2.327903     1.193165      Alps
# 3     0.4826007   0.9897760    -1.212341     2.001142  Caucasus
# 4     0.7651304   1.1635345    -1.622017     1.901734     Chile
# 5     0.1473315   1.3264718    -2.327903     1.973734   Ecuador
# 6    -0.1962311   1.4008344    -1.622017     2.001142 Greenland


coord_radar <- function(...) {
  structure(coord_polar(...), class = c("radar", "polar", "coord"))
}
is.linear.radar <- function(coord) TRUE

##pH
essaiph<-melt(metadata_ph_raw, id.vars = "Site_c",measure.vars = c("pH_median"))
dcastessaiph <-cast(essaiph, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(10, 0), Alps = c(10, 0), Caucasus = c(10, 0),
  Chile = c(10, 0), Ecuador = c(10, 0), Greenland = c(10, 0),
  kyrgyzstan = c(10, 0), Nepal = c(10, 0), New_Zealand = c(10, 0), Norway= c(10, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastessaiph$variable <- NULL
dcastph<- rbind(max_min, dcastessaiph)
rownames(dcastph) <- c("Q3", "Q1","Med")
colnames(dcastph) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                       "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")

create_beautiful_radarchart <- function(data, color = "#00008B", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )

}
##Store plot as object to bind later for final figure
ph_variable<-create_beautiful_radarchart(dcastph, caxislabels = c(0, 2.5, 5, 7.5, 10))
phvariable <- recordPlot()
plot(0) 
replayPlot(phvariable)

##water_temp
watertemp<-melt(metadata_temp_raw, id.vars = "Site_c",measure.vars = c("water_temp_median"))
dcastwatert <-cast(watertemp, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(3.2, 0), Alps = c(3.2, 0), Caucasus = c(3.2, 0),
  Chile = c(3.2, 0), Ecuador = c(3.2, 0), Greenland = c(3.2, 0),
  kyrgyzstan = c(3.2, 0), Nepal = c(3.2, 0), New_Zealand = c(3.2, 0), Norway= c(3.2, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastwatert$variable <- NULL
dcastwt<- rbind(max_min, dcastwatert)
rownames(dcastwt) <- c("Q3", "Q1","Med")
colnames(dcastwt) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                       "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")

op <- par(mar = c(1, 2, 2, 1))
wt_variable<-create_beautiful_radarchart(dcastwt, caxislabels = c(0, 0.8, 1.6, 2.4, 3.2))

##conductivity
condu<-melt(metadata_cond_raw, id.vars = "Site_c",measure.vars = c("conductivity_median"))
dcastcond <-cast(condu, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(104, 0), Alps = c(104, 0), Caucasus = c(104, 0),
  Chile = c(104, 0), Ecuador = c(104, 0), Greenland = c(104, 0),
  kyrgyzstan = c(104, 0), Nepal = c(104, 0), New_Zealand = c(104, 0), Norway= c(104, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastcond$variable <- NULL
dcastcond<- rbind(max_min, dcastcond)
rownames(dcastcond) <- c("Q3", "Q1","Med")
colnames(dcastcond) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                       "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")


##Store plot as object to bind later for final figure
cond_variable<-create_beautiful_radarchart(dcastcond, caxislabels = c(0, 26, 52, 78, 104))
condvariable <- recordPlot()
plot(0)
replayPlot(condvariable)

##gl_sa
glsa<-melt(metadata_gl_sa_raw, id.vars = "Site_c",measure.vars = c("gl_sa_median"))
dcastglsa <-cast(glsa, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(6.4, 0), Alps = c(6.4, 0), Caucasus = c(6.4, 0),
  Chile = c(6.4, 0), Ecuador = c(6.4, 0), Greenland = c(6.4, 0),
  kyrgyzstan = c(6.4, 0), Nepal = c(6.4, 0), New_Zealand = c(6.4, 0), Norway= c(6.4, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastglsa$variable <- NULL
dcastglsa<- rbind(max_min, dcastglsa)
rownames(dcastglsa) <- c("Q3", "Q1","Med")
colnames(dcastglsa) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                         "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")

op <- par(mar = c(1, 2, 2, 1))

##Store plot as object to bind later for final figure
glsa_variable<-create_beautiful_radarchart(dcastglsa, caxislabels = c(0, 1.6, 3.2, 4.8, 6.4))
glsavariable <- recordPlot()
plot(0)
replayPlot(glsavariable)

##turbidity
turb<-melt(metadata_turb_raw, id.vars = "Site_c",measure.vars = c("turb_median"))
dcasturb <-cast(turb, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(760, 0), Alps = c(760, 0), Caucasus = c(760, 0),
  Chile = c(760, 0), Ecuador = c(760, 0), Greenland = c(760, 0),
  kyrgyzstan = c(760, 0), Nepal = c(760, 0), New_Zealand = c(760, 0), Norway= c(760, 0)
)
rownames(max_min) <- c("Max", "Min")
dcasturb$variable <- NULL
dcasturb<- rbind(max_min, dcasturb)
rownames(dcasturb) <- c("Q3", "Q1","Med")
colnames(dcasturb) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                         "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")


##Store plot as object to bind later for final figure
turb_variable<-create_beautiful_radarchart(dcasturb, caxislabels = c(0, 190, 380, 570, 760))
turbvariable <- recordPlot()
plot(0)
replayPlot(turbvariable)


##DOC
doc<-melt(metadata_doc_raw, id.vars = "Site_c",measure.vars = c("doc_median"))
dcastdoc <-cast(doc, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(264, 0), Alps = c(264, 0), Caucasus = c(264, 0),
  Chile = c(264, 0), Ecuador = c(264, 0), Greenland = c(264, 0),
  kyrgyzstan = c(264, 0), Nepal = c(264, 0), New_Zealand = c(264, 0), Norway= c(264, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastdoc$variable <- NULL
dcastdoc<- rbind(max_min, dcastdoc)
rownames(dcastdoc) <- c("Q3", "Q1","Med")
colnames(dcastdoc) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                        "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")


##Store plot as object to bind later for final figure
doc_variable<-create_beautiful_radarchart(dcastdoc, caxislabels = c(0, 66, 132, 198, 264))
docvariable <- recordPlot()
plot(0)
replayPlot(docvariable)

##DIN
din<-melt(metadata_DIN_raw, id.vars = "Site_c",measure.vars = c("DIN_median"))
dcastdin <-cast(din, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(430, 0), Alps = c(430, 0), Caucasus = c(430, 0),
  Chile = c(430, 0), Ecuador = c(430, 0), Greenland = c(430, 0),
  kyrgyzstan = c(430, 0), Nepal = c(430, 0), New_Zealand = c(430, 0), Norway= c(430, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastdin$variable <- NULL
dcastdin<- rbind(max_min, dcastdin)
rownames(dcastdin) <- c("Q3", "Q1","Med")
colnames(dcastdin) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                        "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")
##Store plot as object to bind later for final figure
din_variable<-create_beautiful_radarchart(dcastdin, caxislabels = c(0, 172, 258, 344, 430))
dinvariable <- recordPlot()
plot(0)
replayPlot(dinvariable)

##SRP
srp<-melt(metadata_srp_raw, id.vars = "Site_c",measure.vars = c("srp_median"))
dcastsrp <-cast(srp, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(20, 0), Alps = c(20, 0), Caucasus = c(20, 0),
  Chile = c(20, 0), Ecuador = c(20, 0), Greenland = c(20, 0),
  kyrgyzstan = c(20, 0), Nepal = c(20, 0), New_Zealand = c(20, 0), Norway= c(20, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastsrp$variable <- NULL
dcastsrp<- rbind(max_min, dcastsrp)
rownames(dcastsrp) <- c("Q3", "Q1","Med")
colnames(dcastsrp) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                        "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")

##Store plot as object to bind later for final figure
srp_variable<-create_beautiful_radarchart(dcastsrp, caxislabels = c(0, 5, 10, 15, 20))
srpvariable <- recordPlot()
plot(0)
replayPlot(srpvariable)

##Chla
chla<-melt(metadata_chl_raw, id.vars = "Site_c",measure.vars = c("Chla_median"))
dcastchla <-cast(chla, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(0.16, 0), Alps = c(0.16, 0), Caucasus = c(0.16, 0),
  Chile = c(0.16, 0), Ecuador = c(0.16, 0), Greenland = c(0.16, 0),
  kyrgyzstan = c(0.16, 0), Nepal = c(0.16, 0), New_Zealand = c(0.16, 0), Norway= c(0.16, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastchla$variable <- NULL
dcastchla<- rbind(max_min, dcastchla)
rownames(dcastchla) <- c("Q3", "Q1","Med")
colnames(dcastchla) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                        "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")

##Store plot as object to bind later for final figure
chla_variable<-create_beautiful_radarchart(dcastchla, caxislabels = c(0, 0.4, 0.8, 1.2, 1.6))
chlavariable <- recordPlot()
plot(0)
replayPlot(chlavariable)

##Clays
clays<-melt(metadata_clays_raw, id.vars = "Site_c",measure.vars = c("Clays_median"))
dcastclays <-cast(clays, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(1, 0), Alps = c(1, 0), Caucasus = c(1, 0),
  Chile = c(1, 0), Ecuador = c(1, 0), Greenland = c(1, 0),
  kyrgyzstan = c(1, 0), Nepal = c(1, 0), New_Zealand = c(1, 0), Norway= c(1, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastclays$variable <- NULL
dcastclays<- rbind(max_min, dcastclays)
rownames(dcastclays) <- c("Q3", "Q1","Med")
colnames(dcastclays) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                         "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")

##Store plot as object to bind later for final figure
clays_variable<-create_beautiful_radarchart(dcastclays, caxislabels = c(0, 0.25, 0.5, 0.75, 1))
claysvariable <- recordPlot()
plot(0)
replayPlot(claysvariable)

##Feldspar
feldspar<-melt(metadata_Feldspar_raw, id.vars = "Site_c",measure.vars = c("Feldspar_median"))
dcastfeldspar <-cast(feldspar, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(0.8, 0), Alps = c(0.8, 0), Caucasus = c(0.8, 0),
  Chile = c(0.8, 0), Ecuador = c(0.8, 0), Greenland = c(0.8, 0),
  kyrgyzstan = c(0.8, 0), Nepal = c(0.8, 0), New_Zealand = c(0.8, 0), Norway= c(0.8, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastfeldspar$variable <- NULL
dcastfeldspar<- rbind(max_min, dcastfeldspar)
rownames(dcastfeldspar) <- c("Q3", "Q1","Med")
colnames(dcastfeldspar) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                          "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")


##Store plot as object to bind later for final figure
feldspar_variable<-create_beautiful_radarchart(dcastfeldspar, caxislabels = c(0, 0.2, 0.4, 0.6, 0.8))
feldsparvariable <- recordPlot()
plot(0)
replayPlot(feldsparvariable)

##Qartz
quartz<-melt(metadata_Quartz_raw, id.vars = "Site_c",measure.vars = c("Quartz_median"))
dcastquartz <-cast(quartz, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(0.2, 0), Alps = c(0.2, 0), Caucasus = c(0.2, 0),
  Chile = c(0.2, 0), Ecuador = c(0.2, 0), Greenland = c(0.2, 0),
  kyrgyzstan = c(0.2, 0), Nepal = c(0.2, 0), New_Zealand = c(0.2, 0), Norway= c(0.2, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastquartz$variable <- NULL
dcastquartz<- rbind(max_min, dcastquartz)
rownames(dcastquartz) <- c("Q3", "Q1","Med")
colnames(dcastquartz) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                             "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")

##Store plot as object to bind later for final figure
quartz_variable<-create_beautiful_radarchart(dcastquartz, caxislabels = c(0, 0.05, 0.1, 0.15, 0.2)) 
quartzvariable <- recordPlot() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot(0)
replayPlot(quartzvariable) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

##Calcite
calcite<-melt(metadata_Calcite_raw, id.vars = "Site_c",measure.vars = c("Calcite_median"))
dcastcalcite <-cast(calcite, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(0.048, 0), Alps = c(0.048, 0), Caucasus = c(0.048, 0),
  Chile = c(0.048, 0), Ecuador = c(0.048, 0), Greenland = c(0.048, 0),
  kyrgyzstan = c(0.048, 0), Nepal = c(0.048, 0), New_Zealand = c(0.048, 0), Norway= c(0.048, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastcalcite$variable <- NULL
dcastcalcite<- rbind(max_min, dcastcalcite)
rownames(dcastcalcite) <- c("Q3", "Q1","Med")
colnames(dcastcalcite) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                           "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")

##Store plot as object to bind later for final figure
calcite_variable<-create_beautiful_radarchart(dcastcalcite, caxislabels = c(0, 0.012, 0.024, 0.036, 0.048)) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
calcitevariable <- recordPlot()
plot(0)
replayPlot(calcitevariable)

##glcov
glcov<-melt(metadata_glcov_raw, id.vars = "Site_c",measure.vars = c("gl_cov_median"))
dcastglcov <-cast(glcov, variable ~ Site_c)
max_min <- data.frame(
  Alaska = c(1, 0), Alps = c(1, 0), Caucasus = c(1, 0),
  Chile = c(1, 0), Ecuador = c(1, 0), Greenland = c(1, 0),
  kyrgyzstan = c(1, 0), Nepal = c(1, 0), New_Zealand = c(1, 0), Norway= c(1, 0)
)
rownames(max_min) <- c("Max", "Min")
dcastglcov$variable <- NULL
dcastglcov<- rbind(max_min, dcastglcov)
rownames(dcastglcov) <- c("Q3", "Q1","Med")
colnames(dcastglcov) <- c("Alaska Range","European Alps","Caucasus Mountains", "Chilean Andes", "Ecuadorian Andes",
                            "Southwest Greenland","Pamir & Tien Shan","Himalayas","Southern Alps", "Scandinavian Mountains")

##Store plot as object to bind later for final figure
glcov_variable<-create_beautiful_radarchart(dcastglcov, caxislabels = c(0, 0.2, 0.4, 0.8, 1)) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
glcovvariable <- recordPlot()
plot(0)
replayPlot(glcovvariable) 

library(gridExtra)
library(gridGraphics)
grid.arrange(glcovvariable, glsavariable, nrow = 1)

plot_grid(glcovvariable, NULL, calcitevariable, align = "hv",labels = c("A", "B"), nrow = 1)


plot5 <- plot_grid(glcovariable, NULL, plot2, rel_widths = c(1, -0.5, 1), align = "hv",
                   labels = c("A", "B"), nrow = 1)
