################################################################################
##
## Script for computing taxonomic and functional diversity and dissimilarity 
## between habitats
## 
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)
library(betapart)
library(dplyr)

# loading data
load(here::here("data", "spatial_sp_occ.RData") )
load(here::here("outputs", "sp_3D_coord.RData") )
load(here::here("outputs", "sp_3D_coord_spatial.RData") )

# computing Euclidean distance between species in the 3D space
sp_dist_funct <- dist(sp_3D_coord)

sp_dist_funct <- dist(sp_3D_coord_spatial)

## computing functional diversity based on Hill numbers ####

####### [PS] Same here, don't  know why it's not working for me. I need to use the specific one (in this case, sp_3D_coord_spatial)

mFD :: 

# number of species, functional richness, dispersion and identity (along 3 axes)
spatial_alpha_FDhill <- mFD::alpha.fd.hill (asb_sp_w = spatial_sp_occ,
                    sp_dist = sp_dist_funct,
                    q = c(0,1,2),
                    tau= "mean")

summary(spatial_alpha_FDhill)

spatial_alpha_FDhill$asb_FD_Hill



## computing functional beta-diversity based on Hill numbers ####

# functional dissimilarity
spatial_beta_FDhill <- mFD::beta.fd.hill (asb_sp_w = spatial_sp_occ,
                                           sp_dist = sp_dist_funct,
                                           q = c(0,1, 2),
                                           tau= "mean",
                                           beta_type = "Jaccard")
spatial_beta_FDhill$

# Then use the mFD::dist.to.df function to ease visualizing result

mFD::dist.to.df(list_dist = list("FDBeta_2" = spatial_beta_FDhill$beta_fd_q$q2))


# summary
lapply(spatial_beta_FDhill, summary)
# => beta q0 null because all most distinct pairs of species are in all assemblages
# => beta q=1 still very low, similar composition 

# exploring data to illsutrate indice behaviour
tau_mean<-mean(sp_dist_funct)
tau_mean # 0.39
test_occ <- spatial_sp_occ[1:2,]
sp_test<-names(which(apply(test_occ,2,max)>0))
dist_test<-as.matrix(sp_dist_funct)[sp_test,sp_test]
test_occ[,sp_test]
unique_test<-names(which(apply(test_occ,2,sum)==1))
dist_unique_test<-dist_test[unique_test,unique_test]
summary(as.dist(dist_unique_test))
# => all distances between species unique to different assemblages < mean(dist) 
# so null beta for q=0


# saving ####
save(sp_dist_funct, file=here::here("outputs/", "sp_dist_funct.RData") )
save(spatial_alpha_FDhill, file=here::here("outputs/", "spatial_alpha_FDhill.RData") )
save(spatial_beta_FDhill, file=here::here("outputs/", "spatial_beta_FDhill.RData") )
