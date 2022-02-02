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
load(here::here("outputs", "sp_3D_coord.RData") ) 

## computing taxonomic and functional diversity ####

# number of species, functional richness, dispersion and identity (along 3 axes)
temporal_fd_nokelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = nokelp_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_nokelp <- temporal_fd_nokelp$functional_diversity_indices


temporal_fd_kelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = kelp_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_kelp <- temporal_fd_kelp$functional_diversity_indices



## computing taxonomic and functional beta-diversity ####

# taxonomic dissimilarity = Jaccard index and its components ----
temporal_beta_taxo_nokelp <- betapart::beta.pair(nokelp_sp_occ, index.family = "jaccard")
temporal_beta_taxo_kelp <- betapart::beta.pair(kelp_sp_occ, index.family = "jaccard")

# functional beta no kelp sites ----
# functional dissimilarity = Jaccard-like index and its components
temporal_beta_func_nokelp <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_occ       = nokelp_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

# list of distance matrices with dissimilarity and its turnover
temporal_beta_nokelp <- list (
  taxo_diss = temporal_beta_taxo_nokelp$beta.jac,
  taxo_turn = temporal_beta_taxo_nokelp$beta.jtu,
  func_diss = temporal_beta_func_nokelp$pairasb_fbd_indices$jac_diss,
  func_turn = temporal_beta_func_nokelp$pairasb_fbd_indices$jac_turn
)

# summary
cbind( min=lapply(temporal_beta_nokelp, min), max=lapply(temporal_beta_nokelp, max) )


# functional beta kelp sites ----
temporal_beta_func_kelp <- mFD::beta.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_occ       = kelp_sp_occ,
  check_input      = TRUE,
  beta_family      = c("Jaccard"),
  details_returned = TRUE)

# list of distance matrices with dissimilarity and its turnover
temporal_beta_kelp <- list (
  taxo_diss = temporal_beta_taxo_kelp$beta.jac,
  taxo_turn = temporal_beta_taxo_kelp$beta.jtu,
  func_diss = temporal_beta_func_kelp$pairasb_fbd_indices$jac_diss,
  func_turn = temporal_beta_func_kelp$pairasb_fbd_indices$jac_turn
)

# summary
cbind( min=lapply(temporal_beta_kelp, min), max=lapply(temporal_beta_kelp, max) )


# saving ####

# trait values and trait coding dataframes ----
save(temporal_fd_nokelp, file=here::here("outputs/", "temporal_fd_nokelp.RData") )
save(temporal_alpha_nokelp, file=here::here("outputs/", "temporal_alpha_nokelp.RData") )

save(temporal_fd_kelp, file=here::here("outputs/", "temporal_fd_kelp.RData") )
save(temporal_alpha_kelp, file=here::here("outputs/", "temporal_alpha_kelp.RData") )

save(temporal_beta_nokelp, file=here::here("outputs/", "temporal_beta_nokelp.RData") )
save(temporal_beta_kelp, file=here::here("outputs/", "temporal_beta_kelp.RData") )
