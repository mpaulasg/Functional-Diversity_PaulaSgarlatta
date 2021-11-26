################################################################################
##
## Script for merging biomass from temporal BRUV surveys 
## 
##  Sébastien Villéger
##
## NB: NOT TO BE ON FINAL GITHUB
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


## temporal survey data from sites without kelp ####

# loading raw data from csv----
nokelp_metadata <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_no_kelp.csv"))
head(nokelp_metadata)  

nokelp_replicates_maxN <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_no_kelp.csv"))
head(nokelp_replicates_maxN)  

# merging metadata and data 
nokelp <- left_join(nokelp_metadata, nokelp_replicates_maxN, by=c("Code"="Site") )
head(nokelp)

# computing total biomass per site for each species ----
nokelp_sites_maxN <- nokelp %>% 
  pivot_longer(cols= contains("_"), names_to="species") %>%
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Year, species) %>%
  summarize(total=sum(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = "species", values_from = "total") %>%
  mutate(Code=paste(Site, Year,sep="_"), .before="Site")

  head(nokelp_sites_maxN)
  
nokelp_replicates_maxN %>% select( contains("_") ) %>% sum()
nokelp_sites_maxN %>% select( contains("_") ) %>% sum()
# => same total biomass

# splitting into 2 tables ----
TemporalBRUV_nokelp_metadata <- nokelp_sites_maxN %>% 
  select(Code, Site, Year)

TemporalBRUV_nokelp_species <- nokelp_sites_maxN %>% 
  select( "Code", contains("_") )
  
## saving as csv ----
write.csv(TemporalBRUV_nokelp_metadata, file=here::here("data", "raw_data", "TemporalBRUV_nokelp_metadata.csv"), row.names = FALSE )
write.csv(TemporalBRUV_nokelp_species, file=here::here("data", "raw_data", "TemporalBRUV_nokelp_species.csv"), row.names = FALSE)


## end of no kelp ####
################################################################################

## temporal survey data from sites with kelp ####

# loading raw data from csv----
kelp_metadata <- read.csv(here::here("from_paula", "TemporalBRUV_species_metadata_kelp.csv"))
head(kelp_metadata)  

kelp_replicates_maxN <- read.csv(here::here("from_paula", "TemporalBRUV_species_maxN_kelp.csv"))
head(kelp_replicates_maxN)  

# merging metadata and data 
kelp <- left_join(kelp_metadata, kelp_replicates_maxN, by=c("Code"="Site") )
head(kelp)

# computing total biomass per site for each species ----
kelp_sites_maxN <- kelp %>% 
  pivot_longer(cols= contains("_"), names_to="species") %>%
  mutate(value=as.numeric(value)) %>%
  group_by(Site, Year, species) %>%
  summarize(total=sum(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = "species", values_from = "total") %>%
  mutate(Code=paste(Site, Year,sep="_"), .before="Site")

head(kelp_sites_maxN)

kelp_replicates_maxN %>% select( contains("_") ) %>% sum()
kelp_sites_maxN %>% select( contains("_") ) %>% sum()
# => same total biomass

# splitting into 2 tables ----
TemporalBRUV_kelp_metadata <- kelp_sites_maxN %>% 
  select(Code, Site, Year)

TemporalBRUV_kelp_species <- kelp_sites_maxN %>% 
  select( "Code", contains("_") )

## saving as csv ----
write.csv(TemporalBRUV_kelp_metadata, file=here::here("data", "raw_data", "TemporalBRUV_kelp_metadata.csv"), row.names = FALSE )
write.csv(TemporalBRUV_kelp_species, file=here::here("data", "raw_data", "TemporalBRUV_kelp_species.csv"), row.names= FALSE)


## end of script ####
################################################################################
