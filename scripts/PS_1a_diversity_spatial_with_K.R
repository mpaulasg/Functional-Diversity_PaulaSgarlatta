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
load(here::here("outputs", "sp_3D_coord_spatial_k.RData") ) # [PS] The original one had sp_3D_coord.Rdata


## computing taxonomic and functional diversity ####

# number of species, functional richness, dispersion and identity (along 3 axes)
spatial_fd_k <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord_spatial_k,
  asb_sp_w         = spatial_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

spatial_alpha_k <- spatial_fd_k$functional_diversity_indices
spatial_alpha_k


## computing taxonomic and functional beta-diversity ####

# taxonomic dissimilarity = Jaccard index and its components
spatial_beta_taxo_k <- betapart::beta.pair(spatial_sp_occ, index.family = "jaccard")

# functional dissimilarity = Jaccard-like index and its components
spatial_beta_func_k <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_3D_coord_spatial_k,
  asb_sp_occ       = spatial_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

# list of distance matrices with dissimilarity and its turnover
spatial_beta_k <- list (
  taxo_diss = spatial_beta_taxo_k$beta.jac,
  taxo_turn = spatial_beta_taxo_k$beta.jtu,
  func_diss = spatial_beta_func_k$pairasb_fbd_indices$jac_diss,
  func_turn = spatial_beta_func_k$pairasb_fbd_indices$jac_turn
)

# summary
cbind( min=lapply(spatial_beta_k, min), max=lapply(spatial_beta_k, max) )


# saving ####

# trait values and trait coding dataframes ----
save(spatial_fd_k, file=here::here("outputs/", "spatial_fd_k.RData") )
save(spatial_alpha_k, file=here::here("outputs/", "spatial_alpha_k.RData") )
save(spatial_beta_k, file=here::here("outputs/", "spatial_beta_k.RData") )
