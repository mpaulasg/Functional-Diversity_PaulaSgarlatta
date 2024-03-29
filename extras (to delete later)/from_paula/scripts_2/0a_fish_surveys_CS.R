################################################################################
##
## Script for preparing fish occurrence datasets
## 
## Code by Camille Magneville, Sébastien villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)
library(dplyr)

## temporal survey data ####

# loading raw data from csv----


# metadata and data from sites that use to have kelp and lost it:
kelp_metadata <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_metadata.csv") )
head(kelp_metadata)
unique(kelp_metadata$Site) # 5 sites
unique(kelp_metadata$Year) # 14 years

kelp_sp_maxN <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_species.csv") ) %>%
  column_to_rownames("Code") %>%
  as.matrix()
head(kelp_sp_maxN)

### Joining both datasets

kelp <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_species.csv") )


# summary of surveys and occurrences data ----
kelp_summary <- mFD::asb.sp.summary(asb_sp_w = kelp_sp_maxN)

# retrieve occurrence matrix:
kelp_sp_occ <- kelp_summary$asb_sp_occ


# dimensions
dim(kelp_sp_occ) # 69 assemblages * 101 species

# names of species
kelp_sp <- colnames(kelp_sp_occ) 
length(kelp_sp) # 101 sp


############## => temporal data ready ####

## spatial survey data ####

# metadata of surveys and fish biomass (average across UVC transects) ----
spatial_metadata <- read.csv(here::here("data", "raw_data",  "SpatialUVC_metadata_site.csv"))
head(spatial_metadata)

spatial_sp_biom <- read.csv(here::here("data", "raw_data", "SpatialUVC_species_biomass_site_average.csv")) %>%
  column_to_rownames("Code") %>% 
  as.matrix()

dim(spatial_sp_biom) # 9 assemblages * 53 species

# summary of surveys and occurrences data ----
spatial_summary <- mFD::asb.sp.summary(asb_sp_w = spatial_sp_biom)

# retrieve occurrence matrix:
spatial_sp_occ <- spatial_summary$asb_sp_occ

# species names
spatial_sp <- colnames(spatial_sp_occ)
length(spatial_sp) # 53

## names of species present in at least one dataset ####
sum(spatial_sp %in% kelp_sp)# 35 species shared 

species_allsurveys <- unique( c(kelp_sp ,  spatial_sp) ) 
length(species_allsurveys) # 119 species


## saving dataframes #####
save(kelp_metadata, file=here::here("data", "kelp_metadata.RData") )
save(kelp_sp_occ, file=here::here("data", "kelp_sp_occ.RData") )
save(kelp_summary, file=here::here("data", "kelp_summary.RData") )

save(spatial_metadata, file=here::here("data", "spatial_metadata.RData") )
save(spatial_sp_occ, file=here::here("data", "spatial_sp_occ.RData") )

save(spatial_summary, file=here::here("data", "spatial_summary.RData") )

save(species_allsurveys, file=here::here("data", "species_allsurveys.RData") )

################### end of code ##############################################################

########## For statistics


rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)
library(dplyr)

## temporal survey data ####

# loading raw data from csv----


# metadata and data from sites that use to have kelp and lost it:
kelp_metadata <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_kelp.csv") )
head(kelp_metadata)
unique(kelp_metadata$Site) # 5 sites
unique(kelp_metadata$Year) # 14 years

kelp_sp_maxN <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_kelp.csv") ) %>%
  column_to_rownames("Code") %>%
  as.matrix()
head(kelp_sp_maxN)


# summary of surveys and occurrences data ----
kelp_summary <- mFD::asb.sp.summary(asb_sp_w = kelp_sp_maxN)

# retrieve occurrence matrix:
kelp_sp_occ_transect <- kelp_summary$asb_sp_occ


# dimensions
dim(kelp_sp_occ_transect) # 204 assemblages * 101 species

# names of species
kelp_sp <- colnames(kelp_sp_occ_transect) 
length(kelp_sp) # 101 sp


############## => temporal data ready ####

## spatial survey data ####

# metadata of surveys and fish biomass (average across UVC transects) ----
spatial_metadata <- read.csv(here::here("from_paula",  "SpatialUVC_metadata_transect.csv"))
head(spatial_metadata)

spatial_sp_biom <- read.csv(here::here("from_paula", "SpatialUVC_species_biomass_transect.csv")) %>%
  column_to_rownames("Code") %>% 
  as.matrix()

dim(spatial_sp_biom) # 36 assemblages * 53 species

# summary of surveys and occurrences data ----
spatial_summary_tranect <- mFD::asb.sp.summary(asb_sp_w = spatial_sp_biom)

# retrieve occurrence matrix:
spatial_sp_occ_transect <- spatial_summary_tranect$asb_sp_occ

# species names
spatial_sp_transect <- colnames(spatial_sp_occ_transect)
length(spatial_sp_transect) # 53

## names of species present in at least one dataset ####
sum(spatial_sp_transect %in% kelp_sp)# 35 species shared 

species_allsurveys <- unique( c(kelp_sp ,  spatial_sp_transect) ) 
length(species_allsurveys) # 119 species


## saving dataframes #####

save(kelp_sp_occ_transect, file=here::here("data", "kelp_sp_occ_transect.RData") )

save(spatial_sp_occ_transect, file=here::here("data", "spatial_sp_occ_transect.RData") )

