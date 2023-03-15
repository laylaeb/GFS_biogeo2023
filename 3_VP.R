##Variance partition
library(here)
library(ggplot2)
library(tidyverse)
library(geosphere)   # distm()
library(vegan)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(svd)
library(propack.svd)
# load functions (source Legendre 2016) - not currently used (need to install package AEM)
'PCNM' <-
  function(matdist, thresh=NULL, dbMEM=FALSE, moran=NULL, all=FALSE, include.zero=FALSE, silent=FALSE)
    #
    # Compute the PCNM or dbMEM eigenfunctions corresponding to
    # all eigenvalues (+, 0, -).
    #    In PCNM computation, the diagonal of D = 0.
    #    In dbMEM, the diagonal of D = 4*threshh.
    #    Distance-based MEM are described in Dray et al. 2006.
    #    The name was abbreviated to db-MEM by PPN & PL (subm.)
    # Input file: distance matrix produced by the function "dist".
    # Computation of the threshold requires a function of the library "ape".
    #
    # Original PCNM function: Stephane Dray, November 11, 2004
# The present version: Pierre Legendre, August 2007, January and March 2009
  {
    require(vegan)
    epsilon <- sqrt(.Machine$double.eps)
    a <- system.time({
      if(is.null(moran)) {
        if(dbMEM) { moran=FALSE } else { moran=TRUE }
      }
      single <- FALSE
      if(moran) {
        # cat("The site coordinates were computed from 'matdist'.",'\n')
        pcoa.xy <- pcoa.all(matdist)
        
        if(is.na(pcoa.xy$values[2]) | (pcoa.xy$values[2] < epsilon)) {
          if(!silent) cat("The sites form a straight line on the map.",'\n')
          xy <- pcoa.xy$vectors
          single <- TRUE
        } else {
          xy <- pcoa.xy$vectors[,1:2]
        }
      }
      
      matdist <- as.matrix(matdist)
      n <- nrow(matdist)
      
      # Truncation of distance matrix
      if(is.null(thresh)) {
        spanning <- vegan::spantree(as.dist(matdist))
        threshh <- max(spanning$dist)
        if(!silent) cat("Truncation level =",threshh+0.000001,'\n')
      } else {
        threshh = thresh
        if(!silent) cat("User-provided truncation threshold =",thresh,'\n')
      }
      matdist[matdist > threshh] <- 4*threshh
      
      if(dbMEM==FALSE) { diagonal <- 0 } else { diagonal <- 4*threshh }
      
      mypcnm.all <- pcoa.all(matdist, diagonal=diagonal, all=all, include.zero=include.zero, rn=rownames(matdist))
      
      # Compute Moran's I
      if(moran) {
        require(AEM)
        if(single) {
          nb <- dnearneigh(matrix(c(xy,rep(0,n)),n,2), 0, (threshh + epsilon))
        } else {
          nb <- dnearneigh(xy, 0, (threshh + epsilon))
        }
        fr.to.pcnm2 <- as.matrix(listw2sn(nb2listw(nb))[,1:2])
        weight.dist.coord.mat <- as.matrix(1-(as.dist(matdist)/(4*threshh))^2)
        weight <- weight.dist.coord.mat[fr.to.pcnm2]
        res <- moran.I.multi(mypcnm.all$vectors, link=fr.to.pcnm2, weight=weight)
        Moran <- res$res.mat[,1:2]
        positive <- rep(FALSE,length(mypcnm.all$values))
        positive[which(Moran[,1] > res$expected)] <- TRUE
        Moran <- cbind(as.data.frame(Moran), positive)
        colnames(Moran) <- c("Moran","p.value","Positive")
      }
    })
    a[3] <- sprintf("%2f",a[3])
    if(!silent) cat("Time to compute PCNMs =",a[3]," sec",'\n')
    if(is.null(thresh)) {
      if(moran) {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, Moran_I=Moran, expected_Moran=res$expected, spanning=spanning, thresh=threshh+0.000001)
      } else {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, spanning=spanning, thresh=threshh+0.000001)
      }
    } else {
      if(moran) {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, Moran_I=Moran, expected_Moran=res$expected, thresh=thresh)
      } else {
        res <- list(values=mypcnm.all$values, vectors=mypcnm.all$vectors, thresh=threshh+0.000001)
      }
    }
    res
  }
'pcoa.all' <- function(D, diagonal=0, all=FALSE, include.zero=FALSE, rn=NULL)
  # Principal coordinate decomposition of a square distance matrix D
  # Get the eigenvectors corresponding to all eigenvalues, positive and negative
  # Pierre Legendre, 2005, 2007
  #
  # D : A distance matrix of class 'dist' or 'matrix'.
  # all : If TRUE, the eigenvectors corresponding to all eigenvalues, positive and negative, are shown in the output list.
  # include.zero : If FALSE (default value), the zero eigenvalues as well as their eigenvectors are excluded from the output list.
  # rn : An optional vector of row names, of length n, for the objects.
{
  epsilon <- sqrt(.Machine$double.eps)
  # replace by:     epsilon <- .Machine$double.eps * 10^2
  D <- as.matrix(D)
  n <- nrow(D)
  D <- D + diag(rep(diagonal,n))
  
  # Gower centring, matrix formula
  One <- matrix(1,n,n)
  mat <- diag(n) - One/n
  Dpr2 <- -0.5 * mat %*% (D^2) %*% mat
  trace <- sum(diag(Dpr2))
  
  # Eigenvalue decomposition
  D.eig <- eigen(Dpr2, symmetric=TRUE)
  rel.values <- D.eig$values/trace
  rel.cum <- cumsum(rel.values)
  if(length(rn)!=0) {
    rownames(D.eig$vectors) <- rn
  } else {
    rownames(D.eig$vectors) <- rownames(D)
  }
  
  # Output the results: k eigenvalues and eigenvectors
  if(all) {
    select <- 1:n
    if(!include.zero) {
      exclude <- which(abs(D.eig$values) < epsilon)
      select <- select[-exclude]
    }
    k <- length(select)
    res <- list(values=D.eig$values[select], rel.values=rel.values[select], rel.cum.values=rel.cum[select], vectors=D.eig$vectors[,select], trace=trace)
    # cat("k =",k,"Select =",select,'\n')
    
  } else {
    
    k <- length(which(D.eig$values > epsilon))        
    weight <- sqrt(D.eig$values[1:k])
    if(k == 1) {
      vectors <- D.eig$vectors[,1]*sqrt(D.eig$values[1])
    } else {
      vectors <- D.eig$vectors[,1:k]%*%diag(weight)
    }
    res <- list(values=D.eig$values[1:k], rel.values=rel.values[1:k], rel.cum.values=rel.cum[1:k], vectors=vectors, trace=trace)
  }
  res
}

# load R data
load("comY_2022127.Rdata") # asv_community_trim
load("Env_2022127.Rdata")   # multi_data_na_filter
load("Geo_2022127.Rdata")   # dist_geo
load("geocoordinates_2022127.Rdata")   # metadata_dist_df

##log transformation
bio.b <- log1p(Y.com)

##c'est ici qu'il faut enlever les NA, pas avant
# rename environmental data for ease of use
env.data <- multi_data_na_filter
##We removed feldspar because of multicollinearity
env.names <- c("water_temp", "pH", "conductivity", "turb",  "doc", "srp", "DIN","chla1", "gl_sa", "gl_cov","Calcite","Feldspar","Quartz")
##We removed Clays as there was a signal of multicollinearity
env.mat <- multi_data_na_filter[,env.names]
##Here we only take the environmental variables, we dont keep the mountain ranges
env.mat <- as.data.frame(sapply(env.mat, as.numeric))

# rename geographic coordinates for ease of us
xyz.dat <- metadata_dist_df

# check data structure ## nb of glaciers
with(env.data, table(Site_c))

# Site_c
# Alaska         Alps     Caucasus        Chile      Ecuador    Greenland Kirghizistan        Nepal  New_Zealand       Norway 
# 15           22           17            9            9            6           16           16           18            8 

#### Nested spatial model ####  modified December 7th 2022
## Analyse of microbial spatial variation among and within regions by means of a two-level spatial model
## The among-region component is modeled by a set of dummy variables (N-1 variables with N the number of regions)
## The within-region component is modeled by a set of db-MEM variables for each region
## The db-MEM variables were arranged in blocks corresponding to each region
## within each block, all sites belonging to other regions received the value 0
## See Borcard et al. 2011 (Num ecol with R), Declerck et al. (2011) and function create.MEM.model()

# creating data.frame to store the dbMEMs
var.data.hier <- env.data
library(geosphere)
##create matrix of geographic coordinates
xy.dat = xyz.dat[,1:2]

## ****************************** ##
## SUBSET = Alps                  ##
## ****************************** ##
i = "Alps"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)


## ****************************** ##
## SUBSET = Greenland             ##
## ****************************** ##
i = "Greenland"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region),distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared  
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 
# R2a <- RsquareAdj(rda.region)$adj.r.squared   # R2a = 0.169
# eigenvals(rda.region)[1:2]/sum(eigenvals(rda.region))*RsquareAdj(rda.region)$adj.r.squared/RsquareAdj(rda.region)$r.squared
# axes.region <- scores(rda.region, choices=c(1:2), display="lc", scaling=1)
# pdf(here("Figures", "Figure_Greeland.pdf"), width=12, height=12)
# PCNM.region(axes.region, 1, "Axe1 (16,9%%)")
# dev.off()

## ****************************** ##
## SUBSET = Caucasus              ##
## ****************************** ##
i = "Caucasus"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region),distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared  
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 
# eigenvals(rda.region)[1:2]/sum(eigenvals(rda.region))*RsquareAdj(rda.region)$adj.r.squared/RsquareAdj(rda.region)$r.squared
# axes.region <- scores(rda.region, choices=c(1:2), display="lc", scaling=1)
# pdf(here("Figures", "Figure_Caucasus"), width=12, height=12)
# PCNM.region(axes.region, 1, "Axe1 (5.1%)")
# dev.off()


## ****************************** ##
## SUBSET = Chile                 ##
## ****************************** ##
i = "Chile"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared   # R2a = [1] -0.1062698
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 
# 
# eigenvals(rda.region)[1:2]/sum(eigenvals(rda.region))*RsquareAdj(rda.region)$adj.r.squared/RsquareAdj(rda.region)$r.squared
# axes.region <- scores(rda.region, choices=c(1:2), display="lc", scaling=1)
# pdf(here("Figures", "Figure_Caucasus"), width=12, height=12)
# PCNM.region(axes.region, 1, "Axe1 (5.1%)")
# dev.off()
# 

## ****************************** ##
## SUBSET = Norway                ##
## ****************************** ##
i = "Norway"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
(sel.pos <- which(dbMEM.temp$values > 0))
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared   # R2a = [1] -0.1062698
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 


## ****************************** ##
## SUBSET = Nepal                 ##
## ****************************** ##
i = "Nepal"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared  
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## SUBSET = Kirghizistan          ##
## ****************************** ##
i = "Kirghizistan"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 


## ****************************** ##
## SUBSET = Ecuador               ##
## ****************************** ##
i = "Ecuador"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
sel.pos <- which(dbMEM.temp$values > 0)
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)


# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## SUBSET = New_Zealand           ##
## ****************************** ##
i = "New_Zealand"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
(sel.pos <- which(dbMEM.temp$values > 0))
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call


## ****************************** ##
## SUBSET = Alaska                ##
## ****************************** ##
i = "Alaska"
var.data.temp <- subset(env.data, Site_c == i)
xy.dat.temp <- subset(xy.dat, env.data$Site_c == i)
xyz.dat.temp <- subset(xyz.dat, env.data$Site_c == i)
dist_geo.temp <- distm(xy.dat.temp, fun=distGeo)
bio.b.temp <- subset(bio.b, env.data$Site_c == i)
dbMEM.temp <- pcnm(dist_geo.temp)
(sel.pos <- which(dbMEM.temp$values > 0))
## Eigenfunctions with positive spatial correlation
dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
dbMEM.region <- cbind(data.frame(code_gl=var.data.temp$code_gl), dbMEM.region)
names(dbMEM.region)[2:(length(sel.pos)+1)] <- paste(i, names(dbMEM.region)[2:(length(sel.pos)+1)], sep=".")
var.data.hier <- merge(var.data.hier, dbMEM.region, by="code_gl", sort=F, all.x=T, nomatch=0)

# dbMEM.region <- as.matrix(dbMEM.temp$vectors[,sel.pos])
# rda.region <- dbrda(bio.b.temp~. + Condition(as.matrix(xy.dat.temp)), data=as.data.frame(dbMEM.region), distance="bray")
# R2a <- RsquareAdj(rda.region)$adj.r.squared 
# mod0 <- rda(bio.b.temp ~ 1 + Condition(as.matrix(xy.dat.temp)) , data=as.data.frame(dbMEM.region))
# rda.region.fwd <- ordistep(mod0, scope=formula(rda.region), direction="forward", perm.max=200)
# rda.region.fwd$call
# 


## ****************************** ##
## Editing variables              ##
## ****************************** ##
# Replace all NA in db-MEM by 0
var.data.hier[is.na(var.data.hier)] <- 0
# reorder to match community matrix
var.data.hier$code_gl == env.data$code_gl
var.data.hier2 <- merge(env.data, var.data.hier[,!(colnames(var.data.hier)%in% env.names)], by=c("code_gl","Site_c"), sort=F, all.x=T, nomatch=0)
var.data.hier2$code_gl == env.data$code_gl
# create a region dummy variable (N-1, with N the number of regions)
region.dum <- as.matrix(model.matrix(~-1+var.data.hier2$Site_c))
region.dum <- region.dum[,1:(dim(region.dum)[2]-1)]
# db-MEM only dataset  
region.dbMEM <- var.data.hier2[,16:dim(var.data.hier2)[2]]



#### 3) Models ####
# 3a) variation among regions at the global scale
# Spatial model
rda1.S <-dbrda(bio.b ~. + Condition(as.matrix(xyz.dat)), data=as.data.frame(region.dum), distance="bray")
anova(rda1.S, step=1000)
# Permutation test for dbrda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ `var.data.hier2$Site_cAlaska` + `var.data.hier2$Site_cAlps` + `var.data.hier2$Site_cCaucasus` + `var.data.hier2$Site_cChile` + `var.data.hier2$Site_cEcuador` + `var.data.hier2$Site_cGreenland` + `var.data.hier2$Site_cKirghizistan` + `var.data.hier2$Site_cNepal` + `var.data.hier2$Site_cNew_Zealand` + Condition(as.matrix(xyz.dat)), data = as.data.frame(region.dum), distance = "bray")
# Df SumOfSqs      F Pr(>F)    
# Model      9   10.553 4.4636  0.001 ***
#   Residual 123   32.312                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
R2a <- RsquareAdj(rda1.S)$adj.r.squared #  0.1671687


# Environmental model
rda1.E <- dbrda(bio.b~.+ Condition(as.matrix(xyz.dat)), data=env.mat, distance="bray")
# vif.cca(rda1.E)
# 2.420331                 1.953880                 1.875861                 1.358752                 1.795321                 2.252998                 1.737772 
# doc                      srp                      DIN                    chla1                    gl_sa                   gl_cov                  Calcite 
# 1.286317                 2.117970                 1.795461                 1.165790                 1.258081                 1.363695                 1.941102 
# Feldspar                   Quartz 
# 3.860876                 1.224951  
# library(corrplot)
# library(RColorBrewer)
# library(PerformanceAnalytics)
# chart.Correlation(as.data.frame(env.mat))

eigenvals(rda1.E)[1:2]/sum(eigenvals(rda1.E))*RsquareAdj(rda1.E)$adj.r.squared/RsquareAdj(rda1.E)$r.squared
# dbRDA1     dbRDA2 
# 0.05044934 0.01837816 
anova(rda1.E, step=1000) 
# Permutation test for dbrda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ water_temp + pH + conductivity + turb + doc + srp + DIN + chla1 + gl_sa + gl_cov + Calcite + Feldspar + Quartz + Condition(as.matrix(xyz.dat)), data = env.mat, distance = "bray")
# Df SumOfSqs      F Pr(>F)    
# Model     13    9.320 2.5431  0.001 ***
#   Residual 119   33.545                  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
R2a <- RsquareAdj(rda1.E)$adj.r.squared # R2a =0.1154393
##Compute forward selection
mod0 <- dbrda(bio.b ~ 1 + Condition(as.matrix(xyz.dat)), data=env.mat, distance="bray")
rda1.E.fwd <- ordistep(mod0, scope=formula(rda1.E), direction="forward", perm.max=200)
rda1.E.fwd$call
# 
# dbrda(formula = bio.b ~ Condition(as.matrix(xyz.dat)) + pH + 
#         conductivity + DIN + gl_cov + water_temp + Feldspar + gl_sa + 
#         chla1 + turb, data = env.mat, distance = "bray")

rda2 = dbrda(bio.b ~  Condition(as.matrix(xyz.dat)) + pH + 
               conductivity + DIN +gl_cov + water_temp + Feldspar + gl_sa + 
               chla1 + turb, data = env.mat, distance = "bray")
anova(rda2, by="terms")
# Permutation test for dbrda under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ Condition(as.matrix(xyz.dat)) + pH + conductivity + DIN + gl_cov + water_temp + Feldspar + gl_sa + chla1 + turb, data = env.mat, distance = "bray")
# Df SumOfSqs       F Pr(>F)    
# pH             1    3.046 10.7430  0.001 ***
#   conductivity   1    0.936  3.3003  0.001 ***
#   DIN            1    0.815  2.8723  0.002 ** 
#   gl_cov         1    0.651  2.2945  0.001 ***
#   water_temp     1    0.668  2.3556  0.001 ***
#   Feldspar       1    0.493  1.7388  0.010 ** 
#   gl_sa          1    0.491  1.7312  0.006 ** 
#   chla1          1    0.449  1.5842  0.018 *  
#   turb           1    0.437  1.5403  0.027 *  
#   Residual     123   34.879                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

anova(rda2, by="axis")
# Permutation test for dbrda under reduced model
# Forward tests for axes
# Permutation: free
# Number of permutations: 999
# 
# Model: dbrda(formula = bio.b ~ Condition(as.matrix(xyz.dat)) + pH + conductivity + DIN + gl_cov + water_temp + Feldspar + gl_sa + chla1 + turb, data = env.mat, distance = "bray")
# Df SumOfSqs       F Pr(>F)    
# dbRDA1     1    3.466 12.2230  0.001 ***
#   dbRDA2     1    1.224  4.3161  0.001 ***
#   dbRDA3     1    0.825  2.9085  0.001 ***
#   dbRDA4     1    0.611  2.1545  0.027 *  
#   dbRDA5     1    0.547  1.9293  0.045 *  
#   dbRDA6     1    0.425  1.4998  0.222    
# dbRDA7     1    0.323  1.1402  0.692    
# dbRDA8     1    0.311  1.0953  0.570    
# dbRDA9     1    0.253  0.8934  0.712    
# Residual 123   34.879                   
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

R2a <- RsquareAdj(rda2)$adj.r.squared # [1] 0.110916
eigenvals(rda2)[1:2]/sum(eigenvals(rda2))*RsquareAdj(rda2)$adj.r.squared/RsquareAdj(rda2)$r.squared
# dbRDA1     dbRDA2 
# 0.05626864 0.01986931 

##then subset environmental variables 
env.mat.adj <- subset(env.mat, select = -c(Quartz,srp,Calcite,doc))

vp1 <- varpart(vegdist(bio.b), region.dum, region.dbMEM, env.mat.adj)
vp1;plot(vp1)


# Partition of squared Bray distance in dbRDA 
# 
# Call: varpart(Y = vegdist(bio.b), X = region.dum, region.dbMEM, env.mat.adj)
# 
# Explanatory tables:
#   X1:  region.dum
# X2:  region.dbMEM
# X3:  env.mat.adj 
# 
# No. of explanatory tables: 3 
# Total variation (SS): 50.099 
# No. of observations: 136 
# 
# Partition table:
#   Df R.square Adj.R.square Testable
# [a+d+f+g] = X1         9  0.32091      0.27240     TRUE
# [b+d+e+g] = X2        59  0.39711     -0.07091     TRUE
# [c+e+f+g] = X3         9  0.20428      0.14745     TRUE
# [a+b+d+e+f+g] = X1+X2 68  0.71802      0.43184     TRUE
# [a+c+d+e+f+g] = X1+X3 18  0.40917      0.31827     TRUE
# [b+c+d+e+f+g] = X2+X3 68  0.59206      0.17804     TRUE
# [a+b+c+d+e+f+g] = All 77  0.76527      0.45366     TRUE
# Individual fractions                                   
# [a] = X1 | X2+X3       9               0.27562     TRUE
# [b] = X2 | X1+X3      59               0.13539     TRUE
# [c] = X3 | X1+X2       9               0.02182     TRUE
# [d]                    0              -0.10480    FALSE
# [e]                    0               0.02405    FALSE
# [f]                    0               0.22713    FALSE
# [g]                    0              -0.12555    FALSE
# [h] = Residuals                        0.54634    FALSE
# Controlling 1 table X                                  
# [a+d] = X1 | X3        9               0.17082     TRUE
# [a+f] = X1 | X2        9               0.50275     TRUE
# [b+d] = X2 | X3       59               0.03059     TRUE
# [b+e] = X2 | X1       59               0.15944     TRUE
# [c+e] = X3 | X1        9               0.04587     TRUE
# [c+f] = X3 | X2        9               0.24895     TRUE