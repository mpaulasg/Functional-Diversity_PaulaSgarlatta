################################################################################
##
## Script for computing taxonomic and functional diversity and dissimilarity 
## between habitats using Hill numbers
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
load(here::here("data", "spatial_metadata.RData"))

# computing Euclidean distance between species in the 3D space
sp_dist_funct <- dist(sp_3D_coord)

## computing functional diversity based on Hill numbers ####

#Taxonomic diversity
  
spatial_tax_Hill   <- alpha.fd.hill (asb_sp_w = spatial_sp_occ,
                                          sp_dist = sp_dist_funct,
                                          q = 2, #Recommended in Chao et al 2019
                                          tau= "min") #This gives you the taxonomic profile

TD_spatial <- spatial_tax_Hill$asb_FD_Hill 

# Functional diversity: number of species, functional richness, dispersion and identity (along 3 axes)

spatial_alpha_FDhill <- mFD::alpha.fd.hill (asb_sp_w = spatial_sp_occ,
                    sp_dist = sp_dist_funct,
                    q = c(0, 1, 2) , 
                    tau= "mean") #This gives you functional diversity, as recommended in Chao et al 2019

FD_spatial <- spatial_alpha_FDhill$asb_FD_Hill


## computing tax/functional beta-diversity based on Hill numbers ####

#Taxonomic dissimilarity

spatial_beta_taxhill <-beta.fd.hill (asb_sp_w = spatial_sp_occ,
                                          sp_dist = sp_dist_funct,
                                          q = c(0,1, 2),
                                          tau= "min",
                                          beta_type = "Jaccard")
  
# Then use the mFD::dist.to.df function to ease visualizing result

TD_beta_spatial <- dist.to.df(spatial_beta_taxhill$beta_fd_q) %>% 
  mutate(Habitat1=sub(".*_", "", x1), Habitat2=sub(".*_", "", x2)) %>% 
  select(Habitat1, Habitat2, q0, q1, q2)

# functional dissimilarity
spatial_beta_FDhill <- mFD::beta.fd.hill (asb_sp_w = spatial_sp_occ,
                                           sp_dist = sp_dist_funct,
                                           q = c(0,1, 2),
                                           tau= "mean",
                                           beta_type = "Jaccard")

# Then use the mFD::dist.to.df function to ease visualizing result

FD_beta_spatial <- dist.to.df(spatial_beta_FDhill$beta_fd_q)%>% 
  mutate(Habitat1=sub(".*_", "", x1), Habitat2=sub(".*_", "", x2)) %>% 
  select(Habitat1, Habitat2, q0, q1, q2)

# saving ####
save(sp_dist_funct, file=here::here("outputs/", "sp_dist_funct_spatial_Hill.RData") )
save(TD_spatial, file=here::here("outputs/", "TD_spatial_Hill.RData") )
save(FD_spatial, file=here::here("outputs/", "FD_spatial_Hill.RData") )
save(TD_beta_spatial, file=here::here("outputs/", "TD_beta_spatial_Hill.RData") )
save(FD_beta_spatial, file=here::here("outputs/", "FD_beta_spatial_Hill.RData") )



#From Seb - [PS] Not sure what info is this giving us?

tau_mean<-mean(sp_dist_funct)
tau_mean # 0.39
test_occ <- spatial_sp_occ[1:2,]
sp_test<-names(which(apply(test_occ,2,max)>0))
dist_test<-as.matrix(sp_dist_funct)[sp_test,sp_test]
test_occ[,sp_test]
unique_test<-names(which(apply(test_occ,2,sum)==1))
dist_unique_test<-dist_test[unique_test,unique_test]
summary(as.dist(dist_unique_test))   # all distances between species unique to different assemblages < mean(dist) 
  
# summary
lapply(spatial_beta_FDhill, summary)
# => beta q0 null because all most distinct pairs of species are in all assemblages
# => beta q=1 still very low, similar composition 



