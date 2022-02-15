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

### Joining both datasets

kelp <- read.csv(here::here("data", "raw_data", "TemporalBRUV_kelp_species.csv") )
nokelp <- read.csv(here::here("data", "raw_data", "TemporalBRUV_nokelp_species.csv") )


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

temporal_sp <- unique(c(kelp_sp, nokelp_sp))

length(temporal_sp) # 124 unique species ##[PS] This is not giving unique species, not sure what is exactly giving?

temporal_sp_kelp <- as.data.frame(setdiff(kelp_sp, nokelp_sp)) #[PS] Sp only in kelp

colnames(temporal_sp_kelp) <- "Species"

temporal_sp_nokelp <- as.data.frame(setdiff(nokelp_sp, kelp_sp)) #[PS] Sp only in no kelp

colnames(temporal_sp_nokelp) <- "Species"

temporal_sp_unique <- bind_rows(temporal_sp_kelp, temporal_sp_nokelp)

temporal_sp_unique <- temporal_sp_unique[order(temporal_sp_unique$Species),] 

length(temporal_sp_unique) #41 unique species - THIS IS CORRECT


############## => temporal data ready ####

## spatial survey data ####

# metadata of surveys and fish biomass (average across UVC transects) ----
spatial_metadata <- read.csv(here::here("data", "raw_data",  "SpatialUVC_metadata_site.csv"))
head(spatial_metadata)

spatial_sp_biom <- read.csv(here::here("data", "raw_data", "SpatialUVC_species_biomass_site_average.csv")) %>%
  column_to_rownames("Code") %>% 
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
sum(spatial_sp %in% temporal_sp)# 36 species shared 

species_allsurveys <- unique( c(temporal_sp ,  spatial_sp) ) 
length(species_allsurveys) # 139 species


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


########### For statistics

## temporal survey data ####

# loading raw data from csv----

# metadata and data from sites (and replicates) without kelp
nokelp_metadata_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_no_kelp.csv") )
head(nokelp_metadata_all)
unique(nokelp_metadata_all$Site) # 4 sites
unique(nokelp_metadata_all$Year) # 14 years

nokelp_sp_maxN_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_no_kelp.csv") ) %>%
  column_to_rownames("Site") %>%
  as.matrix()
head(nokelp_sp_maxN_all)

# metadata and data from sites that use to have kelp and lost it:
kelp_metadata_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_kelp.csv") )
head(kelp_metadata_all)
unique(kelp_metadata_all$Site) # 5 sites
unique(kelp_metadata$Year) # 14 years

kelp_sp_maxN_all <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_kelp.csv") ) %>%
  column_to_rownames("Site") %>%
  as.matrix()
head(kelp_sp_maxN_all)


# summary of surveys and occurrences data ----
kelp_summary_all <- mFD::asb.sp.summary(asb_sp_w = kelp_sp_maxN_all)
nokelp_summary_all <- mFD::asb.sp.summary(asb_sp_w = nokelp_sp_maxN_all)

# retrieve occurrence matrix:
kelp_sp_occ_all <- kelp_summary_all$asb_sp_occ
nokelp_sp_occ_all <- nokelp_summary_all$asb_sp_occ


# dimensions
dim(kelp_sp_occ_all) # 204 assemblages * 101 species
dim(nokelp_sp_occ_all) # 167 assemblages * 106 species

############## => temporal data ready ####

## spatial survey data ####

# metadata of surveys and fish biomass (average across UVC transects) ----
spatial_metadata_all <- read.csv(here::here("from_paula",  "SpatialUVC_metadata_transect.csv"))
head(spatial_metadata_all)

spatial_sp_biom_all <- read.csv(here::here("from_paula", "SpatialUVC_species_biomass_transect.csv")) %>%
  column_to_rownames("Site") %>% 
  as.matrix()

dim(spatial_sp_biom_all) # 36 assemblages * 51 species

# summary of surveys and occurrences data ----
spatial_summary_all <- mFD::asb.sp.summary(asb_sp_w = spatial_sp_biom_all)

# retrieve occurrence matrix:
spatial_sp_occ_all <- spatial_summary_all$asb_sp_occ

## saving dataframes #####
save(kelp_metadata_all, file=here::here("data", "kelp_metadata_all.RData") )
save(kelp_sp_occ_all, file=here::here("data", "kelp_sp_occ_all.RData") )
save(kelp_summary_all, file=here::here("data", "kelp_summary_all.RData") )

save(nokelp_metadata_all, file=here::here("data", "nokelp_metadata_all.RData") )
save(nokelp_sp_occ_all, file=here::here("data", "nokelp_sp_occ_all.RData") )
save(nokelp_summary_all, file=here::here("data", "nokelp_summary_all.RData") )

save(spatial_metadata_all, file=here::here("data", "spatial_metadata_all.RData") )
save(spatial_sp_occ_all, file=here::here("data", "spatial_sp_occ_all.RData") )
save(spatial_summary_all, file=here::here("data", "spatial_summary_all.RData") )

############################################### end of script ########################################################
