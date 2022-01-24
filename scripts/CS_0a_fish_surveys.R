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
kelp_nokelp_metadata <- read.csv(here::here("data", "raw_data", "kelp_nokelp_metadata.csv") )

kelp_nokelp <- full_join(kelp, nokelp)

#Changing NA's to 0

kelp_nokelp[is.na(kelp_nokelp)] <- 0

#Adding code to rownames

kelp_nokelp <- kelp_nokelp %>% 
  column_to_rownames("Code") %>% 
  as.matrix()

# summary of surveys and occurrences data ----
kelp_summary <- mFD::asb.sp.summary(asb_sp_w = kelp_sp_maxN)
nokelp_summary <- mFD::asb.sp.summary(asb_sp_w = nokelp_sp_maxN)
kelp_nokelp_summary <- mFD::asb.sp.summary(asb_sp_w = kelp_nokelp)

# retrieve occurrence matrix:
kelp_sp_occ <- kelp_summary$asb_sp_occ
nokelp_sp_occ <- nokelp_summary$asb_sp_occ
kelp_nokelp_occ <- kelp_nokelp_summary$asb_sp_occ


# dimensions
dim(kelp_sp_occ) # 69 assemblages * 101 species
dim(nokelp_sp_occ) # 56 assemblages * 106 species
dim(kelp_nokelp_occ)# 125 assemblages #124 species

# names of species
kelp_sp <- colnames(kelp_sp_occ) #101 sp
length(kelp_sp) # 101 sp

nokelp_sp <- colnames(nokelp_sp_occ) 
length(nokelp_sp) # 106 sp

kelp_nokelp_sp <- colnames(kelp_nokelp_occ) 
length(kelp_nokelp_sp) # 124 sp

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

save(kelp_nokelp_metadata, file=here::here("data", "kelp_nokelp_metadata.RData") )
save(kelp_nokelp_occ, file=here::here("data", "kelp_nokelp_occ.RData") )
save(kelp_nokelp_summary, file=here::here("data", "kelp_nokelp_summary.RData") )

save(spatial_metadata, file=here::here("data", "spatial_metadata.RData") )
save(spatial_sp_occ, file=here::here("data", "spatial_sp_occ.RData") )

save(spatial_summary, file=here::here("data", "spatial_summary.RData") )

save(species_allsurveys, file=here::here("data", "species_allsurveys.RData") )


## end of script ####
