################################################################################
##
## Script for preparing fish occurence datasets
## 
## Code by Camille Magneville, Sébastien villéger and Paula Sgarlatta
##
################################################################################

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)

rm(list=ls()) # cleaning memory

# libraries
library(tidyverse)
library(here)
library(mFD)


## temporal survey data ####

# loading raw data from csv----

# metadata and data from sites without kelp
nokelp_metadata <- read.csv(here::here("data", "raw_data", "TemporalBRUV_nokelp_metadata.csv") )
head(nokelp_metadata)
unique(nokelp_metadata$Site) # 4 sites
unique(nokelp_metadata$Year) # 14 years

nokelp_sp_maxN <- read.csv(here::here("data", "raw_data", "TemporalBRUV_nokelp_species.csv") ) %>%
  column_to_rownames("Code") %>%
  as.matrix()
head(nokelp_sp_maxN)

# metadata and data from sites that use to have kelp and lost it:
kelp_metadata <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_metadata.csv") )
head(kelp_metadata)
unique(kelp_metadata$Site) # 5 sites
unique(kelp_metadata$Year) # 14 years

kelp_sp_maxN <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_species.csv") ) %>%
  column_to_rownames("Code") %>%
  as.matrix()
head(kelp_sp_maxN)


# summary of surveys and occurrences data ----
kelp_summary <- mFD::asb.sp.summary(asb_sp_w = kelp_sp_maxN)
nokelp_summary <- mFD::asb.sp.summary(asb_sp_w = nokelp_sp_maxN)

# retrieve occurrence matrix:
kelp_sp_occ <- kelp_summary$asb_sp_occ
nokelp_sp_occ <- nokelp_summary$asb_sp_occ

# dimensions
dim(kelp_sp_occ) # 69 assemblages * 101 species
dim(nokelp_sp_occ) # 56 assemblages * 106 species

# names of species
kelp_sp <- colnames(kelp_sp_occ)
length(kelp_sp) # 101 sp
nokelp_sp <- colnames(nokelp_sp_occ)
length(nokelp_sp) # 106 sp

sum(kelp_sp %in% nokelp_sp)# 83 species shared

temporal_sp <- unique( c( kelp_sp, nokelp_sp ) )
length(temporal_sp) # 124 unique species

############## => temporal data ready ####

## spatial survey data ####

# metadata of surveys and fish biomass (average across UVC transects) ----
spatial_metadata <- read.csv(here::here("data", "raw_data",  "SpatialUVC_metadata_site.csv"))
head(spatial_metadata)

spatial_sp_biom <- read.csv(here::here("data", "raw_data", "SpatialUVC_species_biomass_site_average.csv"),
                            header = TRUE, row.names = 1) %>%
  as.matrix()
dim(spatial_sp_biom) # 9 assemblages * 51 species

# summary of surveys and occurrences data ----
spatial_summary <- mFD::asb.sp.summary(asb_sp_w = spatial_sp_biom)

# retrieve occurrence matrix:
spatial_sp_occ <- spatial_summary$asb_sp_occ

# species names
spatial_sp <- colnames(spatial_sp_occ)
length(spatial_sp) # 51

## names of species present in at least one dataset ####
sum(spatial_sp %in% temporal_sp)# 32 species shared
species_allsurveys <- unique( c(temporal_sp ,  spatial_sp) )
length(species_allsurveys) # 143 species

## saving dataframes #####
save(kelp_metadata, file=here::here("data", "kelp_metadata.RData") )
save(kelp_sp_occ, file=here::here("data", "kelp_sp_occ.RData") )
save(kelp_summary, file=here::here("data", "kelp_summary.RData") )

save(nokelp_metadata, file=here::here("data", "nokelp_metadata.RData") )
save(nokelp_sp_occ, file=here::here("data", "nokelp_sp_occ.RData") )
save(nokelp_summary, file=here::here("data", "nokelp_summary.RData") )

save(spatial_metadata, file=here::here("data", "spatial_metadata.RData") )
save(spatial_sp_occ, file=here::here("data", "spatial_sp_occ.RData") )
save(spatial_summary, file=here::here("data", "spatial_summary.RData") )

save(species_allsurveys, file=here::here("data", "species_allsurveys.RData") )

## end of script ####
