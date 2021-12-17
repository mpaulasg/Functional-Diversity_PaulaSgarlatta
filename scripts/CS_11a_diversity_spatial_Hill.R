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

# loading data
load(here::here("data", "spatial_sp_occ.RData") )
load(here::here("outputs", "sp_3D_coord.RData") )

# computing Euclidean distance between species in the 3D space
sp_dist_funct <- dist(sp_3D_coord)

## computing functional diversity based on Hill numbers ####

# number of species, functional richness, dispersion and identity (along 3 axes)
spatial_alpha_FDhill <- mFD::alpha.fd.hill (asb_sp_w = spatial_sp_occ,
                    sp_dist = sp_dist_funct,
                    q = c(0,1),
                    tau= "mean",
                    details_returned = FALSE
                    )
summary(spatial_alpha_FDhill)


## computing functional beta-diversity based on Hill numbers ####

# functional dissimilarity
spatial_beta_FDhill <- mFD::beta.fd.hill (asb_sp_w = spatial_sp_occ,
                                           sp_dist = sp_dist_funct,
                                           q = c(0,1),
                                           tau= "mean",
                                           details_returned = FALSE
                                          )

# summary
lapply(spatial_beta_FDhill, summary)
# => beta q0 null because all most distinct pairs of species are in all assemblages
# => beta q=1 still very low, similar composition 

# saving ####
save(sp_dist_funct, file=here::here("outputs/", "sp_dist_funct.RData") )
save(spatial_alpha_FDhill, file=here::here("outputs/", "spatial_alpha_FDhill.RData") )
save(spatial_beta_FDhill, file=here::here("outputs/", "spatial_beta_FDhill.RData") )
