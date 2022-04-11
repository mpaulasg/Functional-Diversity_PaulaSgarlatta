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



# loading data
load(here::here("data", "using biomass-maxN", "kelp_sp_maxN.RData") )
#load(here::here("data", "using biomass-maxN", "kelp_sp_occ.RData") )
load(here::here("outputs", "sp_3D_coord.RData") ) 
load(here::here("data", "using biomass-maxN", "spatial_sp_biom.RData") )
#load(here::here("data", "using biomass-maxN", "spatial_sp_occ.RData") )



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

# saving ####

save(temporal_fd_kelp, file=here::here("outputs/", "using biomass-maxN", "temporal_fd_kelp_biomass.RData") )
save(temporal_alpha_kelp, file=here::here("outputs/","using biomass-maxN",  "temporal_alpha_kelp_biomass.RData") )

save(spatial_fd, file=here::here("outputs/", "using biomass-maxN",  "spatial_fd_biomass.RData") )
save(spatial_alpha, file=here::here("outputs/", "using biomass-maxN", "spatial_alpha_biomass.RData") )


#################################### end of script ####################################################################