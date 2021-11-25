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


## computing taxonomic and functional diversity ####

# number of species, functional richness, dispersion and identity (along 3 axes)
spatial_fd <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = spatial_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

spatial_alpha <- spatial_fd$functional_diversity_indices
spatial_alpha


## computing taxonomic and functional beta-diversity ####

# taxonomic dissimilarity = Jaccard index and its components
spatial_beta_taxo <- betapart::beta.pair(spatial_sp_occ, index.family = "jaccard")

# functional dissimilarity = Jaccard-like index and its components
spatial_beta_func <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_occ       = spatial_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

# list of distance matrices with dissimilarity and its turnover
spatial_beta <- list (
  taxo_diss = spatial_beta_taxo$beta.jac,
  taxo_turn = spatial_beta_taxo$beta.jtu,
  func_diss = spatial_beta_func$pairasb_fbd_indices$jac_diss,
  func_turn = spatial_beta_func$pairasb_fbd_indices$jac_turn
  )

# summary
cbind( min=lapply(spatial_beta, min), max=lapply(spatial_beta, max) )


# saving ####

# trait values and trait coding dataframes ----
save(spatial_fd, file=here::here("outputs/", "spatial_fd.RData") )
save(spatial_alpha, file=here::here("outputs/", "spatial_alpha.RData") )
save(spatial_beta, file=here::here("outputs/", "spatial_beta.RData") )

