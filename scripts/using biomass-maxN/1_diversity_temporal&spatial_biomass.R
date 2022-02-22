################################################################################
##
## Script for computing taxonomic and functional diversity between years for kelp sites,
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
library(betapart)


# loading data
load(here::here("data", "using biomass-maxN", "kelp_sp_maxN.RData") )
load(here::here("data", "using biomass-maxN", "kelp_sp_occ.RData") )
load(here::here("outputs", "sp_3D_coord.RData") ) 
load(here::here("data", "using biomass-maxN", "spatial_sp_biom.RData") )
load(here::here("data", "using biomass-maxN", "spatial_sp_occ.RData") )



## computing taxonomic and functional diversity for kelp sites ####

temporal_fd_kelp <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = kelp_sp_maxN,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_kelp <- temporal_fd_kelp$functional_diversity_indices

## computing taxonomic and functional diversity for spatial data ####

spatial_fd <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = spatial_sp_biom,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

spatial_alpha <- spatial_fd$functional_diversity_indices

## computing taxonomic and functional beta-diversity - spatial ####

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
spatial_beta$taxo_diss
# summary
cbind( min=lapply(spatial_beta, min), max=lapply(spatial_beta, max) )



## computing taxonomic and functional beta-diversity - kelp ####

# taxonomic dissimilarity = Jaccard index and its components ----

temporal_beta_taxo_kelp <- betapart::beta.pair(kelp_sp_occ, index.family = "jaccard")

kelp_turnover <- temporal_beta_taxo_kelp$beta.jtu


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

save(temporal_fd_kelp, file=here::here("outputs/", "using biomass-maxN", "temporal_fd_kelp_biomass.RData") )
save(temporal_alpha_kelp, file=here::here("outputs/","using biomass-maxN",  "temporal_alpha_kelp_biomass.RData") )

save(spatial_fd, file=here::here("outputs/", "using biomass-maxN",  "spatial_fd_biomass.RData") )
save(spatial_alpha, file=here::here("outputs/", "using biomass-maxN", "spatial_alpha_biomass.RData") )

save(temporal_beta_kelp, file=here::here("outputs/", "using biomass-maxN", "temporal_beta_kelp.RData") )
save(spatial_beta, file=here::here("outputs/", "using biomass-maxN", "spatial_beta.RData") )


#################################### end of script ####################################################################

# For statistics

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)


# loading data
load(here::here("data", "using biomass-maxN", "kelp_sp_maxN_transect.RData") )
load(here::here("outputs", "sp_3D_coord.RData") ) 
load(here::here("data", "using biomass-maxN", "spatial_sp_biom_transect.RData") )

## computing taxonomic and functional diversity for kelp sites ####

temporal_fd_kelp_transect <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = kelp_sp_maxN_transect,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

temporal_alpha_kelp_transect <- temporal_fd_kelp_transect$functional_diversity_indices

## computing taxonomic and functional diversity for spatial data ####

spatial_fd_transect <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_3D_coord,
  asb_sp_w         = spatial_sp_biom_transect,
  ind_vect         = c("fide", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

spatial_alpha_transect <- spatial_fd_transect$functional_diversity_indices

# saving ####

# trait values and trait coding dataframes ----

save(temporal_fd_kelp_transect, file=here::here("outputs/","using biomass-maxN", "temporal_fd_kelp_transect.RData") )
save(temporal_alpha_kelp_transect, file=here::here("outputs/","using biomass-maxN", "temporal_alpha_kelp_transect.RData") )

save(spatial_fd_transect, file=here::here("outputs/", "using biomass-maxN", "spatial_fd_transect.RData") )
save(spatial_alpha_transect, file=here::here("outputs/", "using biomass-maxN","spatial_alpha_transect.RData") )
