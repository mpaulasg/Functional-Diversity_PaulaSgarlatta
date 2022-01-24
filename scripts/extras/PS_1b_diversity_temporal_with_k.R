################################################################################
##
## Script for computing taxonomic and functional diversity and dissimilarity 
## between years for kelp and no kelp sites
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
load(here::here("data", "kelp_sp_occ.RData") )
load(here::here("data", "nokelp_sp_occ.RData") )
load(here::here("outputs", "sp_3D_coord_kelp_k.RData") ) # [PS] The original one had sp_3D_coord.Rdata
load(here::here("outputs", "sp_3D_coord_nokelp_k.RData") ) 

## computing taxonomic and functional diversity ####

# number of species, functional richness, dispersion and identity (along 3 axes)
temporal_fd_nokelp_k <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord_nokelp_k,
  asb_sp_w         = nokelp_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_nokelp_k <- temporal_fd_nokelp_k$functional_diversity_indices
temporal_alpha_nokelp_k

temporal_fd_kelp_k <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord_kelp_k,
  asb_sp_w         = kelp_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_kelp_k <- temporal_fd_kelp_k$functional_diversity_indices
temporal_alpha_kelp_k

## computing taxonomic and functional beta-diversity ####

# taxonomic dissimilarity = Jaccard index and its components ----
temporal_beta_taxo_nokelp_k <- betapart::beta.pair(nokelp_sp_occ, index.family = "jaccard")
temporal_beta_taxo_kelp_k <- betapart::beta.pair(kelp_sp_occ, index.family = "jaccard")

# functional beta no kelp sites ----
# functional dissimilarity = Jaccard-like index and its components
temporal_beta_func_nokelp_k <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_3D_coord_nokelp_k,
  asb_sp_occ       = nokelp_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

# list of distance matrices with dissimilarity and its turnover
temporal_beta_nokelp_k <- list (
  taxo_diss_k = temporal_beta_taxo_nokelp_k$beta.jac,
  taxo_turn_k = temporal_beta_taxo_nokelp_k$beta.jtu,
  func_diss_k = temporal_beta_func_nokelp_k$pairasb_fbd_indices$jac_diss,
  func_turn_k = temporal_beta_func_nokelp_k$pairasb_fbd_indices$jac_turn
)

# summary
cbind( min=lapply(temporal_beta_nokelp_k, min), max=lapply(temporal_beta_nokelp_k, max) )


# functional beta kelp sites ----
temporal_beta_func_kelp_k <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_3D_coord_kelp_k,
  asb_sp_occ       = kelp_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

# list of distance matrices with dissimilarity and its turnover
temporal_beta_kelp_k <- list (
  taxo_diss_k = temporal_beta_taxo_kelp_k$beta.jac,
  taxo_turn_k = temporal_beta_taxo_kelp_k$beta.jtu,
  func_diss_k = temporal_beta_func_kelp_k$pairasb_fbd_indices$jac_diss,
  func_turn_k = temporal_beta_func_kelp_k$pairasb_fbd_indices$jac_turn
)

# summary
cbind( min=lapply(temporal_beta_kelp_k, min), max=lapply(temporal_beta_kelp_k, max) )


# saving ####

# trait values and trait coding dataframes ----
save(temporal_fd_nokelp_k, file=here::here("outputs/", "temporal_fd_nokelp_k.RData") )
save(temporal_alpha_nokelp_k, file=here::here("outputs/", "temporal_alpha_nokelp_k.RData") )

save(temporal_fd_kelp_k, file=here::here("outputs/", "temporal_fd_kelp_k.RData") )
save(temporal_alpha_kelp_k, file=here::here("outputs/", "temporal_alpha_kelp_k.RData") )

save(temporal_beta_nokelp_k, file=here::here("outputs/", "temporal_beta_nokelp_k.RData") )
save(temporal_beta_kelp_k, file=here::here("outputs/", "temporal_beta_kelp_k.RData") )
