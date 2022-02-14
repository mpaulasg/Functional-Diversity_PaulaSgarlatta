################################################################################
##
## Script for computing taxonomic and functional diversity between years for kelp and no kelp sites,
##
## and also for spatial data
## 
## Code by Camille Magneville, Sébastien Villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)


# loading data
load(here::here("data", "kelp_sp_occ.RData") )
load(here::here("data", "nokelp_sp_occ.RData") )
load(here::here("outputs", "sp_3D_coord.RData") ) 
load(here::here("data", "spatial_sp_occ.RData") )
load(here::here("outputs", "sp_3D_coord.RData") )

## computing taxonomic and functional diversity for no kelp sites ####

# number of species, functional richness, dispersion and identity (along 3 axes)

temporal_fd_nokelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = nokelp_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_nokelp <- temporal_fd_nokelp$functional_diversity_indices

## computing taxonomic and functional diversity for kelp sites ####

temporal_fd_kelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = kelp_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_kelp <- temporal_fd_kelp$functional_diversity_indices

## computing taxonomic and functional diversity for spatial data ####

spatial_fd <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = spatial_sp_occ,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

spatial_alpha <- spatial_fd$functional_diversity_indices

# saving ####

# trait values and trait coding dataframes ----
save(temporal_fd_nokelp, file=here::here("outputs/", "temporal_fd_nokelp.RData") )
save(temporal_alpha_nokelp, file=here::here("outputs/", "temporal_alpha_nokelp.RData") )

save(temporal_fd_kelp, file=here::here("outputs/", "temporal_fd_kelp.RData") )
save(temporal_alpha_kelp, file=here::here("outputs/", "temporal_alpha_kelp.RData") )

save(spatial_fd, file=here::here("outputs/", "spatial_fd.RData") )
save(spatial_alpha, file=here::here("outputs/", "spatial_alpha.RData") )


#################################### end of script ####################################################################
